#!/bin/bash

############################################################
#        RSV A/B GENOME ASSEMBLY PIPELINE
############################################################

start_total=$(date +%s)

RED='\033[0;31m'
GREEN='\033[1;32m'
YELLOW='\033[1;33m'
CYAN='\033[1;36m'
NC='\033[0m'

############################################################
# HEADER
############################################################
clear
echo -e "${CYAN}"
echo "======================================================"
echo "        RSV A/B GENOME ASSEMBLY PIPELINE"
echo "======================================================"
echo -e "${NC}"

printf "%-15s : %s\n" "Hostname" "$(hostname)"
printf "%-15s : %s cores\n" "CPU" "$(nproc)"
printf "%-15s : %s\n" "RAM" "$(free -h | awk '/Mem:/ {print $2}')"
printf "%-15s : %s\n" "Date" "$(date)"
echo "------------------------------------------------------"

############################################################
# FUNCTIONS
############################################################

progress_bar () {
current=$1
total=$2
label=$3

percent=$((100*current/total))
filled=$((percent/2))
empty=$((50-filled))

bar=$(printf "%${filled}s" | tr ' ' '█')
space=$(printf "%${empty}s")

printf "\r%-25s [%s%s] %3d%% (%d/%d)" \
"$label" "$bar" "$space" \
"$percent" "$current" "$total"
}

step_start(){ STEP_TIME=$(date +%s); }

step_end(){
END=$(date +%s)
DUR=$((END-STEP_TIME))
printf "\n${GREEN}✔ Completed (%02dm %02ds)${NC}\n\n" \
$((DUR/60)) $((DUR%60))
}

############################################################
# CONFIG
############################################################
config=$1
source "$config"

mkdir -p work

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

echo "Samples detected : $total"
echo ""

############################################################
# STEP 1 FASTP
############################################################
echo -e "${YELLOW}[ STEP 1 ] FASTP TRIMMING${NC}"
step_start

c=0
for r1 in "${arr[@]}"; do
((c++))

sample=$(basename "$r1" | sed 's/_R1.*//')
r2=${r1/_R1/_R2}

fastp \
-i "$r1" \
-I "$r2" \
-o work/${sample}_R1.trim.fastq.gz \
-O work/${sample}_R2.trim.fastq.gz \
-h work/${sample}.html \
-j work/${sample}.json \
>/dev/null 2>&1

progress_bar $c $total "Trimming"
done
echo ""
step_end

############################################################
# STEP 2 RSV TYPING
############################################################
echo -e "${YELLOW}[ STEP 2 ] RSV TYPE DETECTION${NC}"
step_start

RSVA=0
RSVB=0
c=0

for r1 in work/*R1.trim.fastq.gz; do
((c++))

sample=$(basename "$r1" _R1.trim.fastq.gz)

bwa mem -t 8 "$RSV_A_REF" \
work/${sample}_R1.trim.fastq.gz \
work/${sample}_R2.trim.fastq.gz \
> work/${sample}_A.sam 2>/dev/null

bwa mem -t 8 "$RSV_B_REF" \
work/${sample}_R1.trim.fastq.gz \
work/${sample}_R2.trim.fastq.gz \
> work/${sample}_B.sam 2>/dev/null

A=$(samtools view -c -F4 work/${sample}_A.sam)
B=$(samtools view -c -F4 work/${sample}_B.sam)

if [ "$A" -gt "$B" ]; then
 mv work/${sample}_A.sam work/${sample}.sam
 rm work/${sample}_B.sam
 echo "$RSV_A_REF" > work/${sample}.ref
 ((RSVA++))
else
 mv work/${sample}_B.sam work/${sample}.sam
 rm work/${sample}_A.sam
 echo "$RSV_B_REF" > work/${sample}.ref
 ((RSVB++))
fi

progress_bar $c $total "Typing RSV"
done
echo ""
step_end

echo -e "${GREEN}RSV Typing Summary${NC}"
echo "RSV-A : $RSVA"
echo "RSV-B : $RSVB"
echo ""

############################################################
# STEP 3 IVAR TRIM
############################################################
if [ "$bedfile" == "Yes" ]; then
echo -e "${YELLOW}[ STEP 3 ] PRIMER TRIMMING${NC}"
step_start

for sam in work/*.sam; do
base=${sam%.sam}
ivar trim -i "$sam" -b "$bedfile_path" -p "$base" \
>/dev/null 2>&1
done

step_end
fi

############################################################
# STEP 4 BAM PROCESS
############################################################
echo -e "${YELLOW}[ STEP 4 ] BAM PROCESSING${NC}"
step_start

for sam in work/*.sam; do
base=${sam%.sam}

samtools view -bS "$sam" |
samtools sort -o ${base}.sorted.bam

samtools index ${base}.sorted.bam
done

step_end

############################################################
# STEP 5 COVERAGE
############################################################
echo -e "${YELLOW}[ STEP 5 ] COVERAGE ANALYSIS${NC}"
step_start

for bam in work/*.sorted.bam; do
samtools coverage "$bam" > ${bam%.bam}_coverage.tsv
samtools depth "$bam" > ${bam%.bam}_depth.tsv
done

step_end

############################################################
# STEP 6 CONSENSUS
############################################################
echo -e "${YELLOW}[ STEP 6 ] CONSENSUS GENERATION${NC}"
step_start

for bam in work/*.sorted.bam; do
base=${bam%.sorted.bam}
ref=$(cat ${base}.ref)

samtools mpileup -aa -A -d 0 -Q 0 \
-f "$ref" "$bam" |
ivar consensus -t 0.6 \
-p ${base}.fa >/dev/null 2>&1
done

step_end

############################################################
# FINAL OUTPUT
############################################################
DATE=$(date +%F)
mkdir -p $DATE

mv work/*.fa $DATE/
mv work/*coverage.tsv $DATE/
mv work/*depth.tsv $DATE/

cat $DATE/*.fa > $DATE/combined.fasta

############################################################
# SUMMARY
############################################################
end_total=$(date +%s)
runtime=$((end_total-start_total))

echo ""
echo "======================================================"
echo -e "${GREEN} PIPELINE COMPLETED SUCCESSFULLY ${NC}"
echo "------------------------------------------------------"
echo "Samples     : $total"
echo "RSV-A       : $RSVA"
echo "RSV-B       : $RSVB"
echo "Output      : $DATE/"
printf "Runtime     : %02dh %02dm %02ds\n" \
$((runtime/3600)) $((runtime%3600/60)) $((runtime%60))
echo "======================================================"
echo ""
