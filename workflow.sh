#!/bin/bash

# Menyimpan waktu mulai skrip
start_time=$(date +%s)

# List Warna
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[1;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Progress Bar
progress_bar() {
    local current_step=$1
    local total_steps=$2
    local step_name=$3
    local progress=$((current_step * 100 / total_steps))
    local bar_length=50
    local filled_length=$((progress * bar_length / 100))
    local empty_length=$((bar_length - filled_length))

    local bar=$(printf '%*s' "$filled_length" | tr ' ' '*')$(printf '%*s' "$empty_length" | tr ' ' '.')
    # Menampilkan progres dalam format (1/16)
    printf "\r\033[1;36m%-40s (%d/%d) [%s] %d%%\033[0m" "$step_name" "$current_step" "$total_steps" "$bar" "$progress"
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
bedfile_enabled="No"
bedfile_path=""

while IFS='=' read -r key value; do
    case "$key" in
        fastq_folder) fastq_folder="$value" ;;
        reference_folder) reference_folder="$value" ;;
        bedfile) bedfile_enabled="$value" ;;  # Membaca konfigurasi bedfile
        bedfile_path) bedfile_path="$value" ;;  # Membaca path bedfile jika ada
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

# Pastikan bedfile ada jika diperlukan
if [ "$bedfile_enabled" == "Yes" ] && [ ! -f "$bedfile_path" ]; then
    echo -e "${RED}Error: Bedfile not found at $bedfile_path.${NC}"
    exit 1
fi
echo "" # Pindah ke baris baru setelah progress selesai
echo "" # Pindah ke baris baru setelah progress selesai
echo "" # Pindah ke baris baru setelah progress selesai
echo -e "${YELLOW}_____________________________________________________________________${NC}"
echo "" # Pindah ke baris baru setelah progress selesai
# Resource Komputer
printf "${CYAN}%-20s${NC}  ${BLUE}%-50s${NC}\n" "CPU Model" "$(lscpu | grep 'Model name' | sed 's/Model name:/ /' | awk '{$1=$1};1')"
printf "${CYAN}%-20s${NC}  ${BLUE}%-50s${NC}\n" "Total RAM" "$(free -h | grep Mem | awk '{print $2}')"
printf "${CYAN}%-20s${NC}  ${BLUE}%-50s${NC}\n" "Sistem Operasi" "$(uname -s) $(uname -r)"
printf "${CYAN}%-20s${NC}  ${BLUE}%-50s${NC}\n" "Arsitektur" "$(uname -m)"
printf "${CYAN}%-20s${NC}  ${BLUE}%-50s${NC}\n" "Jumlah Core CPU" "$(nproc) cores"
printf "${CYAN}%-20s${NC}  ${BLUE}%-50s${NC}\n" "Hostname" "$(hostname)"
echo -e "${YELLOW}_____________________________________________________________________${NC}"
echo "" # Pindah ke baris baru setelah progress selesai
echo "" # Pindah ke baris baru setelah progress selesai

# 3. Fastp Trimming Step
total_files=$(find "$fastq_folder"/*.fastq.gz | wc -l)
arr=( $(ls "$fastq_folder"/*.fastq.gz) )
total_samples=$((total_files / 2))
current_sample=0

echo -e "${YELLOW}Step 1: Quality Control (Fastp) ${NC}"

for ((i=0; i<$total_files; i+=2)); do
    ((current_sample++))
    sample_name=$(basename "${arr[$i]}" .fastq.gz)

    # Run fastp for trimming
    (fastp -i "${arr[$i]}" -I "${arr[$i+1]}" -o "${sample_name}_trimmed_R1.fastq.gz" -O "${sample_name}_trimmed_R2.fastq.gz" -h "${sample_name}_fastp.html" -j "${sample_name}_fastp.json" 2>/dev/null) &
    pid=$!
    
    while ps -p $pid > /dev/null; do
        progress_bar $current_sample $total_samples "Trimming with fastp"
        sleep 0.5
    done
done

echo -e "${GREEN}\nCompleted!${NC}"

# 4. Proses Mapping
echo "" # Pindah ke baris baru setelah progress selesai
echo -e "${YELLOW}Step 2: Read to Reference (BWA) ${NC}"

# Updated to make sure progress is calculated correctly
total_samples=$((total_files / 2))
current_sample=0

for ((i=0; i<$total_files; i+=2)); do
    ((current_sample++))
    sample_name=$(basename "${arr[$i]}" .fastq.gz)

    # Run BWA MEM with trimmed reads
    (bwa mem -t 12 "$ref_file" "${sample_name}_trimmed_R1.fastq.gz" "${sample_name}_trimmed_R2.fastq.gz" > "$sample_name.sam" 2>/dev/null) &
    pid=$!
    
    while ps -p $pid > /dev/null; do
        progress_bar $current_sample $total_samples "Mapping to Reference"
        sleep 0.5
    done
done

echo -e "${GREEN}\nCompleted!${NC}"

# Continue with other steps like trimming with ivar, SAM to BAM conversion, sorting, etc.
echo "" # Pindah ke baris baru setelah progress selesai

# 4. ivar trim (Langkah Terpisah)
if [ "$bedfile_enabled" == "Yes" ]; then
    echo -e "${YELLOW}Step 3: Trimming Primer (iVar) ${NC}"

    sam_files=$(ls *.sam)
    total_trim=$(echo "$sam_files" | wc -l)
    current_trim=0

    for sam_file in $sam_files; do
        ((current_trim++))
        progress_bar $current_trim $total_trim "Trimming Primer (ivar) "

        # Langkah ivar trim
        (ivar trim -i "$sam_file" -b "$bedfile_path" -p "${sam_file%.sam}" 2>/dev/null ) &
        pid=$!
        
        while ps -p $pid > /dev/null; do
            progress_bar $current_trim $total_trim "Trimming with ivar"
            sleep 0.5
        done
    done

    echo -e "${GREEN}\nCompleted!${NC}\n"
fi

# 5. Konversi SAM ke BAM
mkdir -p "$fastq_folder/temporary"
mv *.sam "$fastq_folder/temporary/"
sam_files=$(ls "$fastq_folder/temporary"/*.sam 2>/dev/null)
total_conversions=$(echo "$sam_files" | wc -l)
current_conversion=0

echo -e "${YELLOW}Step 4: Convert SAM to BAM (Samtools) ${NC}"

for sam_file in $sam_files; do
    ((current_conversion++))
    progress_bar $current_conversion $total_conversions "Convert SAM to BAM"
    
    # Menggunakan samtools view untuk konversi SAM ke BAM
    (samtools view -h -F 4 -F 2048 -@ 7 "$sam_file" > "${sam_file%.sam}.bam" 2>/dev/null) &
    pid=$!

    while ps -p $pid > /dev/null; do
        progress_bar $current_conversion $total_conversions "Convert SAM to BAM"
        sleep 0.5
    done
done
echo -e "${GREEN}\nCompleted!${NC}\n"

# 6. Sorting BAM files
bam_files=$(ls "$fastq_folder/temporary"/*.bam 2>/dev/null)
total_sorts=$(echo "$bam_files" | wc -l)
current_sort=0

echo -e "${YELLOW}Step 5: Sorting BAM files (Samtools) ${NC}"

for bam_file in $bam_files; do
    ((current_sort++))
    progress_bar $current_sort $total_sorts "Sorting BAM files"
    
    (samtools sort -@ 7 "$bam_file" > "${bam_file%.bam}_sorted.bam" 2>/dev/null) &
    pid=$!

    while ps -p $pid > /dev/null; do
        progress_bar $current_sort $total_sorts "Sorting BAM files"
        sleep 0.5
    done
done

echo -e "${GREEN}\nCompleted!${NC}\n"

# 7. Menghitung Coverage dan Depth
sorted_bam_files=$(ls "$fastq_folder/temporary"/*_sorted.bam 2>/dev/null)
total_coverages=$(echo "$sorted_bam_files" | wc -l)
current_coverage=0

echo -e "${YELLOW}Step 6: Coverage Analysis (Samtools) ${NC}"

for sorted_bam_file in $sorted_bam_files; do
    ((current_coverage++))
    progress_bar $current_coverage $total_coverages "Count Coverage"
    
    (samtools coverage "$sorted_bam_file" > "${sorted_bam_file%.bam}_coverage.tsv" 2>/dev/null) &
    pid=$!

    while ps -p $pid > /dev/null; do
        progress_bar $current_coverage $total_coverages "Count Coverage"
        sleep 0.5
    done
done

echo -e "${GREEN}\nCompleted!${NC}\n"

# Menghitung Depth
total_depths=$(echo "$sorted_bam_files" | wc -l)
current_depth=0

echo -e "${YELLOW}Step 7: Depth Analysis (Samtools) ${NC}"

for sorted_bam_file in $sorted_bam_files; do
    ((current_depth++))
    progress_bar $current_depth $total_depths "Count Depth"
    
    (samtools depth "$sorted_bam_file" > "${sorted_bam_file%.bam}_depth.tsv" 2>/dev/null) &
    pid=$!

    while ps -p $pid > /dev/null; do
        progress_bar $current_depth $total_depths "Count Depth"
        sleep 0.5
    done
done

echo -e "${GREEN}\nCompleted!${NC}"

echo "" # Pindah ke baris baru setelah progress selesai
# 8. Membuat file konsensus
total_consensus=$(echo "$sorted_bam_files" | wc -l)
current_consensus=0

echo -e "${YELLOW}Step 8: Generating Consensus(iVar) ${NC}"

for sorted_bam_file in $sorted_bam_files; do
    ((current_consensus++))
    progress_bar $current_consensus $total_consensus "Generating Consensus"
    
    (samtools mpileup -aa -A -d 0 -Q 0 -q 0 "$sorted_bam_file" 2>/dev/null | ivar consensus -t 0.2 -p  "${sorted_bam_file%.bam}.fa" > /dev/null 2>&1) &
    pid=$!

    while ps -p $pid > /dev/null; do
        progress_bar $current_consensus $total_consensus "Generating Consensus"
        sleep 0.5
    done
done

echo -e "${GREEN}\nCompleted!${NC}"

# 9. Output folder berdasarkan tanggal
current_date=$(date +"%Y-%m-%d")
mkdir -p "$current_date"
mv "$fastq_folder/temporary"/*_coverage.tsv "$current_date/"
mv "$fastq_folder/temporary"/*_depth.tsv "$current_date/"
mv "$fastq_folder/temporary"/*.fa "$current_date/"
rm *.fastq.gz
rm *.bam
rm *.html*
rm *.json
cat "$current_date"/*.fa > "$current_date"/combined.fasta

# Menghitung waktu selesai skrip
end_time=$(date +%s)

# Menghitung durasi dalam detik
duration=$((end_time - start_time))

# Mengonversi detik ke format jam:menit:detik
hours=$((duration / 3600))
minutes=$(((duration % 3600) / 60))
seconds=$((duration % 60))

echo ""
echo -e "${GREEN}Workflow completed in ${hours} hours, ${minutes} minutes, and ${seconds} seconds.${NC}"
echo ""
