#!/bin/bash

# Warna untuk output yang lebih menarik
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[1;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Fungsi untuk menampilkan progress bar dengan animasi
progress_bar() {
    local current_step=$1
    local total_steps=$2
    local step_name=$3
    local progress=$((current_step * 100 / total_steps))
    local bar_length=50
    local filled_length=$((progress * bar_length / 100))
    local empty_length=$((bar_length - filled_length))

    local bar=$(printf '%*s' "$filled_length" | tr ' ' '*')$(printf '%*s' "$empty_length" | tr ' ' '.')
    printf "\r\033[1;36m%-40s %d/%d [%s] %d%%\033[0m" "$step_name" "$current_step" "$total_steps" "$bar" "$progress"
}

# 1. Baca file konfigurasi
if [ -z "$1" ]; then
    echo -e "${RED}Error: No config file provided. Please provide a config file.${NC}"
    exit 1
fi

config_file=$1
if [ ! -f "$config_file" ]; then
    echo -e "${RED}Error: Config file not found.${NC}"
    exit 1
fi

# Set variables dari file konfigurasi
while IFS='=' read -r key value; do
    case "$key" in
        fastq_folder) fastq_folder="$value" ;;
        reference_folder) reference_folder="$value" ;;
    esac
done < "$config_file"

# 2. Validasi folder referensi dan FASTQ
if [ ! -d "$reference_folder" ]; then
    echo -e "${RED}Error: Reference folder not found.${NC}"
    exit 1
fi

ref_file=$(find "$reference_folder" -name "*.fasta" | head -n 1)
if [ -z "$ref_file" ]; then
    echo -e "${RED}Error: No .fasta file found in the reference folder.${NC}"
    exit 1
fi

if [ ! -d "$fastq_folder" ]; then
    echo -e "${RED}Error: FASTQ folder not found.${NC}"
    exit 1
fi

echo -e "${YELLOW}Pipeline ini dibuat oleh Markus Evan Anggia

       /\       
      /  \      /\     
     /    \    /  \    
    /      \  /    \   /\    
   /        \/      \ /  \   
  /                  /    \  
 /__________________/      \ 
  |   ________  __  |  _  _| 
  |  |  _     ||  | | | || | 
  |  | | |    ||  | | | || | 
  |__|_|_|____||__|_|_|_||_|_|


${NC}"

echo -e "${YELLOW}${NC}"



# 3. Proses Mapping
total_files=$(find "$fastq_folder"/*.fastq.gz | wc -l)
arr=( $(ls "$fastq_folder"/*.fastq.gz) )
total_samples=$((total_files / 2))
current_sample=0

echo -e "${YELLOW}Step 1: Mapping to Reference using BWA MEM${NC}"

for ((i=0; i<$total_files; i+=2)); do
    ((current_sample++))
    sample_name=$(basename "${arr[$i]}" .fastq.gz)

    (bwa mem -t 12 "$ref_file" "${arr[$i]}" "${arr[$i+1]}" > "$sample_name.sam" 2>/dev/null) &
    pid=$!
    
    while ps -p $pid > /dev/null; do
        progress_bar $current_sample $total_samples "Mapping to Reference"
        sleep 0.5
    done
    echo "" # Pindah ke baris baru setelah progress selesai
done

echo -e "${GREEN}\nMapping completed!${NC}\n"

# 4. Konversi SAM ke BAM
mkdir -p "$fastq_folder/temporary"
mv *.sam "$fastq_folder/temporary/"
sam_files=$(ls "$fastq_folder/temporary"/*.sam 2>/dev/null)
total_conversions=$(echo "$sam_files" | wc -l)
current_conversion=0

echo -e "${YELLOW}Step 2: Convert SAM to BAM${NC}"

for sam_file in $sam_files; do
    ((current_conversion++))
    progress_bar $current_conversion $total_conversions "Convert SAM to BAM"
    
    (samtools view -h -F 4 -F 2048 -@ 7 "$sam_file" > "${sam_file%.sam}.bam" 2>/dev/null) &
    pid=$!

    while ps -p $pid > /dev/null; do
        progress_bar $current_conversion $total_conversions "Convert SAM to BAM"
        sleep 0.5
    done
    echo "" # Pindah ke baris baru setelah progress selesai
done
echo -e "${GREEN}\nConversion completed!${NC}\n"

# 5. Sorting BAM files
bam_files=$(ls "$fastq_folder/temporary"/*.bam 2>/dev/null)
total_sorts=$(echo "$bam_files" | wc -l)
current_sort=0

echo -e "${YELLOW}Step 3: Sorting BAM files${NC}"

for bam_file in $bam_files; do
    ((current_sort++))
    progress_bar $current_sort $total_sorts "Sorting BAM files"
    
    (samtools sort -@ 7 "$bam_file" > "${bam_file%.bam}_sorted.bam" 2>/dev/null) &
    pid=$!

    while ps -p $pid > /dev/null; do
        progress_bar $current_sort $total_sorts "Sorting BAM files"
        sleep 0.5
    done
    echo "" # Pindah ke baris baru setelah progress selesai
done

echo -e "${GREEN}\nSorting completed!${NC}\n"

# 6. Menghitung Coverage dan Depth
sorted_bam_files=$(ls "$fastq_folder/temporary"/*_sorted.bam 2>/dev/null)
total_coverages=$(echo "$sorted_bam_files" | wc -l)
current_coverage=0

echo -e "${YELLOW}Step 4: Coverage Analysis${NC}"

for sorted_bam_file in $sorted_bam_files; do
    ((current_coverage++))
    progress_bar $current_coverage $total_coverages "Count Coverage"
    
    (samtools coverage "$sorted_bam_file" > "${sorted_bam_file%.bam}_coverage.tsv" 2>/dev/null) &
    pid=$!

    while ps -p $pid > /dev/null; do
        progress_bar $current_coverage $total_coverages "Count Coverage"
        sleep 0.5
    done
    echo "" # Pindah ke baris baru setelah progress selesai
done

echo -e "${GREEN}\nCoverage Analysis completed!${NC}\n"

# Menghitung Depth
total_depths=$(echo "$sorted_bam_files" | wc -l)
current_depth=0

echo -e "${YELLOW}Step 5: Depth Analysis${NC}"

for sorted_bam_file in $sorted_bam_files; do
    ((current_depth++))
    progress_bar $current_depth $total_depths "Count Depth"
    
    (samtools depth "$sorted_bam_file" > "${sorted_bam_file%.bam}_depth.tsv" 2>/dev/null) &
    pid=$!

    while ps -p $pid > /dev/null; do
        progress_bar $current_depth $total_depths "Count Depth"
        sleep 0.5
    done
    echo "" # Pindah ke baris baru setelah progress selesai
done

echo -e "${GREEN}\nDepth Analysis completed!${NC}\n"

# 7. Membuat file konsensus
total_consensus=$(echo "$sorted_bam_files" | wc -l)
current_consensus=0

echo -e "${YELLOW}Step 6: Generating Consensus${NC}"

for sorted_bam_file in $sorted_bam_files; do
    ((current_consensus++))
    progress_bar $current_consensus $total_consensus "Generating Consensus"
    
    (samtools mpileup -aa -A -d 0 -Q 0 -q 0 "$sorted_bam_file" 2>/dev/null | ivar consensus -p "${sorted_bam_file%.bam}.fasta" > /dev/null 2>&1) &
    pid=$!

    while ps -p $pid > /dev/null; do
        progress_bar $current_consensus $total_consensus "Generating Consensus"
        sleep 0.5
    done
    echo "" # Pindah ke baris baru setelah progress selesai
done

echo -e "${GREEN}\nConsensus generation completed!${NC}"

# 8. Output folder berdasarkan tanggal
current_date=$(date +"%Y-%m-%d")
mkdir -p "$current_date"
mv "$fastq_folder/temporary"/*_coverage.tsv "$current_date/"
mv "$fastq_folder/temporary"/*_depth.tsv "$current_date/"
mv "$fastq_folder/temporary"/*.fa "$current_date/"

echo -e "${GREEN}\nAll steps completed! Output saved in folder: $current_date${NC}"
