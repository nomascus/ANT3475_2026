#!/bin/bash  
#SBATCH --job-name=bwa_pipeline  
#SBATCH --cpus-per-task=1  
#SBATCH --mem=1G  
#SBATCH --time=00:05:00  
#SBATCH --output=%x_%j.out  
#SBATCH --error=%x_%j.err

############################
# 0. Load modules
############################

module load bwa
module load samtools
module load bcftools

############################
# 1. Define input files.
############################

REF=$1
R1=$2
R2=$3
SAMPLE=$4

############################
# 2. Index reference files # don't need to run
############################

# echo "Indexing reference..."

# samtools faidx $REF 
# bwa index $REF

############################
# 3. Align reads to reference
############################

echo "Aligning reads..."

bwa mem $REF $R1 $R2 > ${SAMPLE}.sam

############################
# 4. Convert SAM → sorted BAM
############################

echo "Converting to BAM..."

samtools view -bS ${SAMPLE}.sam | \
samtools sort -o ${SAMPLE}.sorted.bam

############################
# 5. Index BAM
############################

samtools index ${SAMPLE}.sorted.bam

############################
# 6. Basic BAM filtering
# (remove low mapping quality reads)
############################

echo "Filtering BAM..."

samtools view -b -q 20 ${SAMPLE}.sorted.bam > ${SAMPLE}.filtered.bam

samtools index ${SAMPLE}.filtered.bam

############################
# 7. Call variants
############################

echo "Calling variants..."

bcftools mpileup -Ou -f $REF ${SAMPLE}.filtered.bam | \
bcftools call -mv -Ov -o ${SAMPLE}.vcf

############################
# 8. Filter VCF
# (keep QUAL > 30)
############################

echo "Filtering VCF..."

bcftools filter -i 'QUAL>30' ${SAMPLE}.vcf -o ${SAMPLE}.filtered.vcf

echo "Pipeline complete!"
