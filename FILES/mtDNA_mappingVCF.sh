#!/bin/bash

############################
# 0. Load modules
############################

module load bwa
module load samtools
module load bcftools

############################
# 1. Define input files
############################

REF=Human_mtDNA.fasta
R1=demo_mt_R1.fastq
R2=demo_mt_R2.fastq
SAMPLE=demo_sample

############################
# 2. Align reads to reference
############################

echo "Aligning reads..."

bwa mem -t 4 $REF $R1 $R2 > ${SAMPLE}.sam

############################
# 3. Convert SAM → sorted BAM
############################

echo "Converting to BAM..."

samtools view -bS ${SAMPLE}.sam | \
samtools sort -o ${SAMPLE}.sorted.bam

############################
# 4. Index BAM
############################

samtools index ${SAMPLE}.sorted.bam

############################
# 5. Basic BAM filtering
# (remove low mapping quality reads)
############################

echo "Filtering BAM..."

samtools view -b -q 20 ${SAMPLE}.sorted.bam > ${SAMPLE}.filtered.bam

samtools index ${SAMPLE}.filtered.bam

############################
# 6. Call variants
############################

echo "Calling variants..."

bcftools mpileup -Ou -f $REF ${SAMPLE}.filtered.bam | \
bcftools call -mv -Ov -o ${SAMPLE}.vcf

############################
# 7. Filter VCF
# (keep QUAL > 30)
############################

echo "Filtering VCF..."

bcftools filter -i 'QUAL>30' ${SAMPLE}.vcf -o ${SAMPLE}.filtered.vcf

echo "Pipeline complete!"
