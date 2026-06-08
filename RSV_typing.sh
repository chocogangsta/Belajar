#!/bin/bash

# List Warna
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[1;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

echo ""
echo ""
echo ""

############################################################
# RSV A/B GENOME ASSEMBLY PIPELINE
############################################################

start_total=$(date +%s)

GREEN='\033[1;32m'
NC='\033[0m'

mkdir -p work logs

############################################################
# CONFIG (TIDAK DIUBAH)
############################################################
config=$1
source "$config"

############################################################
# SYSTEM INFO
############################################################
HOSTNAME=$(hostname)
CPU=$(nproc)
RAM=$(free -h | awk '/Mem:/ {print $2}')
DATE_NOW=$(date)

############################################################
# HEADER
############################################################
echo "======================================================"
echo "        RSV A/B GENOME ASSEMBLY PIPELINE"
echo "======================================================"
echo ""
printf "%-15s : %s\n" "Hostname" "$HOSTNAME"
printf "%-15s : %s cores\n" "CPU" "$CPU"
printf "%-15s : %s\n" "RAM" "$RAM"
printf "%-15s : %s\n" "Date" "$DATE_NOW"
echo "------------------------------------------------------"

############################################################
# PROGRESS BAR
############################################################
progress_bar(){

cur=$1
tot=$2
label=$3

percent=$((100*cur/tot))
fill=$((percent/2))
empty=$((50-fill))

bar=$(printf "%${fill}s" | tr ' ' '#')
space=$(printf "%${empty}s")

printf "\r%-25s [%s%s] %3d%% (%d/%d)" \
"$label" "$bar" "$space" \
"$percent" "$cur" "$tot"
}

############################################################
# SUMMARY FILE
############################################################
TYPING_SUMMARY="work/rsv_typing.txt"
COVERAGE_SUMMARY="work/coverage.txt"
> $TYPING_SUMMARY
> $COVERAGE_SUMMARY

############################################################
# REFERENCES
############################################################
RSV_A_REF=$(find "$reference_folder" -name "*A*.fasta")
RSV_B_REF=$(find "$reference_folder" -name "*B*.fasta")

bwa index "$RSV_A_REF" >/dev/null 2>&1
bwa index "$RSV_B_REF" >/dev/null 2>&1

############################################################
# FASTQ
############################################################
arr=( $(ls $fastq_folder/*R1*.fastq.gz) )
total=${#arr[@]}
echo ""
echo "Samples detected : $total"
echo ""

############################################################
# STEP 1 FASTP
############################################################
echo -e "${YELLOW}[ STEP 1 ] FASTP TRIMMING${NC}"
step_start=$(date +%s)

i=0
for r1 in "${arr[@]}"; do
((i++))

sample=$(basename "$r1" | sed 's/_R1.*//')
r2=${r1/_R1/_R2}

fastp -q 20 -A -i "$r1" -I "$r2" \
-o work/${sample}_R1.trim.fastq.gz \
-O work/${sample}_R2.trim.fastq.gz \
>/dev/null 2>&1

progress_bar $i $total "Trimming"
done
echo ""

step_end=$(date +%s)
printf "\n${GREEN}✔ Completed (%02dm %02ds)${NC}\n\n" \
$(((step_end-step_start)/60)) \
$(((step_end-step_start)%60))

############################################################
# STEP 2 RSV TYPE DETECTION
############################################################
echo -e "${YELLOW}[ STEP 2 ] RSV TYPE DETECTION${NC}"

RSVA=0
RSVB=0
step_start=$(date +%s)

i=0
for r1 in work/*R1.trim.fastq.gz; do
((i++))

sample=$(basename "$r1" _R1.trim.fastq.gz)

bwa mem "$RSV_A_REF" work/${sample}_R1.trim.fastq.gz \
work/${sample}_R2.trim.fastq.gz \
> work/${sample}_A.sam 2>/dev/null

bwa mem "$RSV_B_REF" work/${sample}_R1.trim.fastq.gz \
work/${sample}_R2.trim.fastq.gz \
> work/${sample}_B.sam 2>/dev/null

A=$(samtools view -q 20 -c -F4 work/${sample}_A.sam)
B=$(samtools view -q 20 -c -F4 work/${sample}_B.sam)

if [ "$A" -gt "$B" ]; then
 mv work/${sample}_A.sam work/${sample}.sam
 rm work/${sample}_B.sam
 echo "$RSV_A_REF" > work/${sample}.ref
 echo "${sample} : RSV-A" >> $TYPING_SUMMARY
 ((RSVA++))
else
 mv work/${sample}_B.sam work/${sample}.sam
 rm work/${sample}_A.sam
 echo "$RSV_B_REF" > work/${sample}.ref
 echo "${sample} : RSV-B" >> $TYPING_SUMMARY
 ((RSVB++))
fi

progress_bar $i $total "Typing RSV"
done
echo ""

step_end=$(date +%s)
printf "\n${GREEN}✔ Completed (%02dm %02ds)${NC}\n\n" \
$(((step_end-step_start)/60)) \
$(((step_end-step_start)%60))

echo "RSV Typing Summary"
cat $TYPING_SUMMARY
echo ""

############################################################
# STEP 4 BAM PROCESSING
############################################################
echo -e "${YELLOW}[ STEP 3 ] BAM PROCESSING${NC}"

step_start=$(date +%s)
i=0

for sam in work/*.sam; do
((i++))

sample=$(basename "$sam" .sam)

samtools view -q 20 -bS "$sam" |
samtools sort -o work/${sample}.sorted.bam >/dev/null 2>&1

samtools index work/${sample}.sorted.bam >/dev/null 2>&1

progress_bar $i $total "Processing BAM"
done
echo ""

step_end=$(date +%s)
printf "\n${GREEN}✔ Completed (%02dm %02ds)${NC}\n\n" \
$(((step_end-step_start)/60)) \
$(((step_end-step_start)%60))

############################################################
# STEP 5 COVERAGE
############################################################
echo -e "${YELLOW}[ STEP 4 ] COVERAGE ANALYSIS${NC}"

step_start=$(date +%s)

for bam in work/*.sorted.bam; do

sample=$(basename "$bam" .sorted.bam)
ref=$(cat work/${sample}.ref)

GENOME_LEN=$(grep -v ">" "$ref" | wc -c)

COV_BASE=$(samtools depth "$bam" \
| awk '$3>=10{c++} END{print c+0}')

COV_PERCENT=$(awk -v a=$COV_BASE -v b=$GENOME_LEN \
'BEGIN{printf "%.0f",(a/b)*100}')

echo "${sample} : ${COV_PERCENT}%" >> $COVERAGE_SUMMARY

done

cat $COVERAGE_SUMMARY
echo ""

step_end=$(date +%s)
printf "${GREEN}✔ Completed (%02dm %02ds)${NC}\n\n" \
$(((step_end-step_start)/60)) \
$(((step_end-step_start)%60))

############################################################
# STEP 6 CONSENSUS
############################################################
echo -e "${YELLOW}[ STEP 5 ] CONSENSUS GENERATION${NC}"

step_start=$(date +%s)
i=0

for bam in work/*.sorted.bam; do
((i++))

sample=$(basename "$bam" .sorted.bam)
ref=$(cat work/${sample}.ref)

samtools mpileup -aa -A -d 0 -Q 0 \
-f "$ref" "$bam"  2>> logs/bcftools.log |
ivar consensus -t 0.6 -q 20 -m 20 \
-p work/${sample}.fa \
>/dev/null 2>&1

progress_bar $i $total "Consensus"
done
echo ""

step_end=$(date +%s)
printf "\n${GREEN}✔ Completed (%02dm %02ds)${NC}\n\n" \
$(((step_end-step_start)/60)) \
$(((step_end-step_start)%60))

############################################################
# OUTPUT
############################################################
DATE=$(date +%F)
mkdir -p $DATE
mv work/*.fa $DATE/
cat $DATE/*.fa > $DATE/combined.fasta

############################################################
# FINAL SUMMARY
############################################################
end_total=$(date +%s)
runtime=$((end_total-start_total))

echo "======================================================"
echo " PIPELINE COMPLETED SUCCESSFULLY"
echo "------------------------------------------------------"
echo "Samples     : $total"
echo "RSV-A       : $RSVA"
echo "RSV-B       : $RSVB"
echo "Output      : $DATE/"
printf "Runtime     : %02dh %02dm %02ds\n" \
$((runtime/3600)) \
$((runtime%3600/60)) \
$((runtime%60))
echo "======================================================"
