# Devoir 5


Download the following files for today's exercises

```
wget https://raw.githubusercontent.com/nomascus/ANT3475_2026/refs/heads/main/FILES/Mapping_VCF/reference_mt.fasta
wget https://raw.githubusercontent.com/nomascus/ANT3475_2026/refs/heads/main/FILES/Mapping_VCF/mt_multiSNP_2000_R1.fastq 
wget https://raw.githubusercontent.com/nomascus/ANT3475_2026/refs/heads/main/FILES/Mapping_VCF/mt_multiSNP_2000_R2.fastq
```

## File Descriptions

| file name | file type |file description
| - | - | - | - |
| reference_mt.fasta | FASTA | human mitochondrial reference genome| 
| mt_multiSNP_2000_R1.fastq | FASTQ | 2000 human mitochondrial reads (Forward reads) |
| mt_multiSNP_2000_R1.fastq | FASTQ | 2000 human mitochondrial reads (Reverse reads)|
 

1) Load the modules for bwa, samtools, and bcftools

2) Make all the index files for the reference genome

3) Map the reads to the reference genome

4) Compress, sort and index the SAM file

5) Filter the BAM File to remove reads with mapping scores less than 20

6) Call variants (make a VCF) from the filtered BAM file

7) Filter the VCF to remove SNPs with quality scores less than 30

8) Which SNP was removed from the VCF during the filtration process in step 7?
