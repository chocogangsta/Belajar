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

# 3. Running TB-Profiler
total_files=$(find "$fastq_folder"/*.fastq.gz | wc -l)
arr=( $(ls "$fastq_folder"/*.fastq.gz) )
total_samples=$((total_files / 2))
current_sample=0

echo -e "${YELLOW}Step 1: Running TB-Profiler ${NC}"

for ((i=0; i<$total_files; i+=2)); do
    ((current_sample++))
    sample_name=$(basename "${arr[$i]}" .fastq.gz)

    # Run fastp for trimming
    (
    bash -c "
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate TBprofiler 
    echo 'Environment aktif: '$CONDA_DEFAULT_ENV >/dev/null 2>&1
    tb-profiler profile -1 '${arr[$i]}' -2 '${arr[$i+1]}' -p '${sample_name}' --txt >/dev/null 2>&1
    conda deactivate
    "
) &
    pid=$!
    


    while ps -p $pid > /dev/null; do
        progress_bar $current_sample $total_samples "Run TB Profiler"
        sleep 0.5
    done
done

echo -e "${GREEN}\nCompleted!${NC}"


# 4. Make Consensus Genome

vcf_files=$(ls vcf/*.targets.vcf.gz 2>/dev/null)
total_sorts=$(echo "$vcf_files" | wc -w)
current_sort=0

for vcf_file in $vcf_files; do
    ((current_sort++))
    progress_bar $current_sort $total_sorts "Make Consensus Genome"

    sample_name=$(basename "$vcf_file" .targets.vcf.gz)

    (
    bcftools index "$vcf_file"
    bcftools consensus -f "$ref_file" -o "vcf/${sample_name}.consensus.fasta" "$vcf_file" 2>/dev/null
    ) &
    pid=$!

    while ps -p $pid > /dev/null; do
        progress_bar $current_sort $total_sorts "Make Consensus Genome"
        sleep 0.5
    done
done

echo -e "${GREEN}\nCompleted!${NC}\n"

echo -e "${CYAN}Step 4: Running Prokka annotation...${NC}"

# Aktifkan environment Prokka
source $(conda info --base)/etc/profile.d/conda.sh
conda activate assembly_env

# Jalankan Prokka untuk setiap hasil TB-Profiler (FASTA hasil assembly)
for fasta in vcf/*.fasta; do
    sample=$(basename "$fasta" .fasta)
    echo -e "${BLUE}Running Prokka for ${sample}...${NC}"

    prokka "$fasta" \
        --outdir prokka_results/${sample} \
        --prefix ${sample} \
        --gcode 11 \
        --proteins genbank/sequence.gb \
        --cpus 4 \
        --force
done


# 9. Output folder berdasarkan tanggal
current_date=$(date +"%Y-%m-%d")
mkdir -p "$current_date"

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
