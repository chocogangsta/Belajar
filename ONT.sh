#!/bin/bash

# Warna
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; CYAN='\033[0;36m'; NC='\033[0m'

# Progress bar
progress_bar() {
    local current=$1; local total=$2; local title=$3
    local percent=$((current * 100 / total))
    local filled=$((percent / 2)); local empty=$((50 - filled))
    local bar=$(printf '%*s' "$filled" | tr ' ' '#')$(printf '%*s' "$empty" | tr ' ' '.')
    printf "\r${CYAN}%-35s (%d/%d) [%s] %d%%%s" "$title" "$current" "$total" "$bar" "$percent" "$NC"
}

# Baca konfigurasi
read_config() {
    config_file=$1
    if [ ! -f "$config_file" ]; then echo -e "${RED}Config file not found!${NC}"; exit 1; fi
    while IFS='=' read -r key value; do
        case "$key" in
            fastq_folder) fastq_folder="$value" ;;
            reference_folder) reference_folder="$value" ;;
            bedfile) bedfile_enabled="$value" ;;
            bedfile_path) bedfile_path="$value" ;;
        esac
    done < "$config_file"
}

# Validasi input
validate_input() {
    [ ! -d "$fastq_folder" ] && echo -e "${RED}FASTQ folder not found!${NC}" && exit 1
    [ ! -d "$reference_folder" ] && echo -e "${RED}Reference folder not found!${NC}" && exit 1
    ref_file=$(find "$reference_folder" -name "*.fasta" | head -n 1)
    [ -z "$ref_file" ] && echo -e "${RED}No .fasta in reference folder!${NC}" && exit 1
    [ "$bedfile_enabled" == "Yes" ] && [ ! -f "$bedfile_path" ] && echo -e "${RED}BED file missing!${NC}" && exit 1
    mkdir -p "$fastq_folder/temporary"
}

run_fastp() {
    echo -e "${YELLOW}Step 1: Trimming with fastp${NC}"
    files=("$fastq_folder"/*.fastq.gz); total=${#files[@]}; count=0
    for file in "${files[@]}"; do
        ((count++)); name=$(basename "$file" .fastq.gz)
        progress_bar $count $total "Trimming: $name"
        fastp -i "$file" -A -q 5 -o "$fastq_folder/temporary/${name}_trimmed.fastq.gz" -h "$fastq_folder/temporary/${name}_fastp.html" >/dev/null 2>&1
    done
    echo -e "${GREEN}\nCompleted fastp${NC}"
}

run_mapping() {
    echo -e "${YELLOW}Step 2: Mapping with minimap2${NC}"
    files=("$fastq_folder/temporary"/*_trimmed.fastq.gz); total=${#files[@]}; count=0
    for file in "${files[@]}"; do
        ((count++)); name=$(basename "$file" _trimmed.fastq.gz)
        progress_bar $count $total "Mapping: $name"
        minimap2 -ax map-ont --secondary=no -t 4 "$ref_file" "$file" > "$fastq_folder/temporary/${name}.sam" 2>/dev/null
    done
    echo -e "${GREEN}\nCompleted mapping${NC}"
}

run_variant_calling() {
    echo -e "${YELLOW}Step 3: Variant Calling + Filtering${NC}"
    for sam in "$fastq_folder/temporary/"*.sam; do
        name=$(basename "$sam" .sam)
        bam="$fastq_folder/temporary/${name}.bam"
        sorted="$fastq_folder/temporary/${name}_sorted.bam"

        samtools view -bS "$sam" > "$bam"
        samtools sort "$bam" -o "$sorted"
        samtools index "$sorted"

        bcftools mpileup -Ou -f "$ref_file" "$sorted" | \
        bcftools call -mv -Ou -V indels | \
        bcftools filter -i 'QUAL>=5 && DP>=5' -Ob -o "$fastq_folder/temporary/${name}_variants.bcf"

        bcftools index "$fastq_folder/temporary/${name}_variants.bcf"
    done
    echo -e "${GREEN}Completed variant filtering${NC}"
}

extract_bcf_to_bed() {
    echo -e "${YELLOW}Step 4: Extracting BCF to BED${NC}"
    for bcf in "$fastq_folder/temporary/"*_variants.bcf; do
        name=$(basename "$bcf" _variants.bcf)
        bcftools query -f '%CHROM\t%POS0\t%END\n' "$bcf" > "$fastq_folder/temporary/${name}_ext.bed"
    done
    echo -e "${GREEN}Completed extracting BED${NC}"
}

create_masking() {
    echo -e "${YELLOW}Step 5: Create Masking BED${NC}"
    for sorted_bam in "$fastq_folder/temporary/"*_sorted.bam; do
        name=$(basename "$sorted_bam" _sorted.bam)
        bed="$fastq_folder/temporary/${name}_ext.bed"
        mask="$fastq_folder/temporary/${name}_mask.bed"
        if [ -f "$bed" ]; then
            bedtools genomecov -bga -ibam "$sorted_bam" | awk '$4 < 5' | \
            bedtools subtract -a - -b "$bed" > "$mask"
        fi
    done
    echo -e "${GREEN}Completed masking${NC}"
}

generate_consensus() {
    echo -e "${YELLOW}Step 6: Generate Consensus Fasta${NC}"
    for mask in "$fastq_folder/temporary/"*_mask.bed; do
        name=$(basename "$mask" _mask.bed)
        bcf="$fastq_folder/temporary/${name}_variants.bcf"
        [ -f "$bcf" ] && bcftools consensus -f "$ref_file" -m "$mask" --haplotype A "$bcf" > "$fastq_folder/temporary/${name}_consensus.fasta"
    done
    echo -e "${GREEN}Consensus fasta generated${NC}"
}

move_output() {
    output_date=$(date +"%Y-%m-%d")
    mkdir -p "$output_date/consensus" "$output_date/report"
    mv "$fastq_folder/temporary/"*_consensus.fasta "$output_date/consensus" 2>/dev/null
    mv "$fastq_folder/temporary/"*fastp.html "$output_date/report" 2>/dev/null
    echo -e "${GREEN}Output moved to: ${output_date}${NC}"
}

main() {
    read_config "$1"
    validate_input
    run_fastp
    run_mapping
    run_variant_calling
    extract_bcf_to_bed
    create_masking
    generate_consensus
    move_output
    echo -e "${GREEN}All steps completed successfully.${NC}"
}

# Jalankan
if [ $# -lt 1 ]; then
    echo -e "${RED}Usage: $0 <config.txt>${NC}"; exit 1
else
    main "$1"
fi
