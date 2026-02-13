
# Mapping reads and calling variants

Files to download for today's class

```
wget https://raw.githubusercontent.com/nomascus/ANT3475_2026/refs/heads/main/FILES/Human_mtDNA.fasta
wget https://raw.githubusercontent.com/nomascus/ANT3475_2026/refs/heads/main/FILES/demo_mt_R1.fastq
wget https://raw.githubusercontent.com/nomascus/ANT3475_2026/refs/heads/main/FILES/demo_mt_R2.fastq
wget https://raw.githubusercontent.com/nomascus/ANT3475_2026/refs/heads/main/FILES/mtDNA_mappingVCF.sh
```

Imagine if you had just downloaded DNA sequencing data from 200 human samples and your goal is to identify how the populations are subdivided by geography. Unfortunately, we can't just start working on the biological questions, because all the DNA sequencing machine gives us is a giant mixed up pile of **short DNA fragments, called "reads"** that come in a file type called FASTQ (note that this different from FASTA). In order to analyze the data, we need to generate a VCF file, which contains only the nucleotide positions that are variable in our sample. in other words need to go from the raw sequencing reads to the SNPs. That's the goal for today.

### What Are FASTQ Files?

FASTQ files contain four lines of text for each sequencing read:

Example:

```
@read1_Forward
ACTGATCGTACGATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
+
IIIIIIIIIIHHHHHHHHHHFFFFFFFFFFDDDDDDDDDD!!!!!!!!!!
@read2_Forward
GCTAGCTAGCATCGATCGATCGTAGCTAGCATCGATCGTAGCTAGCATCG
+
IIIIIIHHHHHHFFFFFFDDDDDD!!!!!!!!!IIIIHHHHFFFFDDDDD

```
This is the raw data from the sequencing machine. It contains:

1) A header starting with some information about the read we are looking at. FASTQ headers always start with the @ character.

2) A line of nucleotide sequence

3) a plus sign

4) A line of characters corresponding to the quality score for each nucleotide. In other words, how confident is the sequencing machine that the first nucleotide in the read is A, the second nucleotide is C, etc. This quality score lets us remove untrustworthy reads from our analysis

*Note*: DNA sequencing is usually done in the forward and the reverse directions along the DNA strand, so we usually get two FASTQ files for each sample, a forward and a reverse file.

Before we can start studying the biological questions, we need to:

1.  Figure out where the reads belong in the genome
2.  Clean up bad data  
3.  Find differences from the reference genome for each individual

# 1. Mapping (FASTQ → SAM)

Before we do anything, we need to load three programs that are stored internally on the server. To do these, we use the **module load** command. If these programs don't run later, check that you have loaded them first!

```
module load samtools
module load bwa
module load bcftools
```

The mapping process takes a FASTQ file and assigns each read to the most likely position in a reference genome. With enough reads, we can cover the entire genome. 

Imagine that you had a shredded copy of a document. It's kind of like taking each piece of shredded paper and lining it up against another copy of the document to find where it belongs. 
 
![Hands-on: Mapping / Mapping / Sequence analysis](https://training.galaxyproject.org/training-material/topics/sequence-analysis/images/mapping/mapping.png)

The difference here is that we use areference genome, which is just a FASTA file of all the nucleotides lined up in all the chromosomes for the organism we want to study. For example, the human genome is 3 billion nucleotides long, and we have to find where each sequencing read fits. 

Because a reference genome file can be huge, we need to build an index file that helps other programs find the mapped reads more quickly. This is analogous to having an index in the back of a cookbook helping you find the page of a recipe.

```
samtools faidx reference.fasta 
bwa index reference.fasta 
```
Here, we use the **samtools index** and **bwa index** commands to make the index files. The index files will automatically be created and end in .amb, .ann, .bwt, .fai, .pac, and .sa  

To map reads, we use a program called **bwa mem**

```
bwa mem -t 1 reference.fasta reads_sample1_F.fastq reads_sample1_R.fastq > mappedReads_sample1.sam
``` 

We call the program bwa mem, then give it three files:

1) The reference genome
2) The forward reads (often compressed in .gz format)
3) The reverse reads (often compressed in .gz format)

and then save it (redirect the output of bwa mem) into a new file, which is called a sam file.

The **-t 1** option indicates that we want to use a single processor on the computer. For a really big file, we might need to increase this number, but this is a tiny example.

### What is a SAM file?

SAM (Sequence Alignment/Map) is a text file that stores:

-   Where each read maps in the genome 
-   The mapping quality
-   Alignment details
  
   Here is an example of what a SAM file looks like. 

  ```
  @HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000

read1	0	chr5	1540	60	50M	*	0	0	ACTGATCGTACGATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII!!!!!!!!!#####
read2	0	chr1	292	15	50M	*	0	0	GCTAGCTAGCATCGATCGATCGTAGCTAGCATCGATCGTAGCTAGCATCG	#####!!!!!!!!!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

  ```
We don't need to go into details, but basically it says: 

Read one should start at position 1540 on chromosome 5, and it is confident that it belongs there (mapping quality score of 60), 

Read two should start at position 292 on chromosome 1, and it is not very confident that it belongs there (mapping quality score of 15), 

In a real SAM file there will be millions of reads all of which cover each nucleotide position many times.

# 2.  Compress, sort and index the SAM file

Now that we have a our mapping file, we need to organize it and compress it, because many downstream programs require reads to be ordered by chromosome and nucleotide position. The sam file we have is huge and the reads are in a random order.  To do this, we use a group of programs called **samtools**.

```
samtools view -bS mappedReads_sample1.sam | samtools sort -o mappedReads_sample1.sorted.bam
```
Let's break this down 
1)  **samtools view** opens the sam file. **-bS** tells it that the input file will be a SAM and the output will be a compressed version called a **BAM** file.

2) We pipe the output into **samtools sort**, which sorts the bam file by chromosome and nucleotide position

3)  **-o** tells us that the name of the output file will be mappedReads_sample1.sorted.bam

There is one more small preparatory step. Because the BAM file can be huge, we need to build an index file that helps other programs find the mapped reads more quickly. 
```
samtools index mappedReads_sample1.sorted.bam
```
Here, we use the **samtools index** command to make the index file. The index file will automatically be created and end in .bai

```
mappedReads_sample1.sorted.bam.bai 
```

# 3. Filter the BAM File

Remember that there is a mapping quality score for each read in our BAM file. The computer does it's best to find a position in the genome for every read, but sometimes it mistakenly maps a read to the wrong place or multiple places in the genome. This can be the result of many factors, such as contaminant bacterial DNA in the sample, degraded DNA, poorly made reference genomes, or sequencing errors. We need to remove the reads that map poorly so we only work with nucleotide sequences that are from the right organism and mapped to the right place.

```
samtools view -b -q 20 mappedReads_sample1.sorted.bam > mappedReads_sample1.sorted.filtered.bam
``` 

Here again, we are using **samtools view** to open  bam file with the -b option, then we add **-q 20**,  which tells removes reads with mapping quality scores below 20. This step greatly improves data reliability. Then we save it as a new file using the redirect symbol **>** 

# 4. Variant Calling (BAM → VCF)

Now that we have a nice clean bam file, it's time to make the **VCF**. Remember that this is the **Variant Call File**,  which gives us a list of all the variable positions in the genome, and the genotype of each individual. in our sample. 

### Calling variants with bcftools 

```
bcftools mpileup -Ou -f reference.fasta mappedReads_sample1.sorted.filtered.bam | \
bcftools call -mv -Ov -o HumanSamples.vcf
``` 
**bcftools** is another suite of programs for working with genomic data. We are going to use two of its programs

1.  `mpileup` looks at every position in the genome and sees how many reads cover that position. For example:

![](https://github.com/nomascus/ANT3475_2026/blob/main/FILES/mpileup.png?raw=true)

**-Ou** means that the output type is an uncompressed BAM

**-f** specifies the reference genome file

Note that here we are just using a single BAM from a single sample. However, if we want to want to feed the program bam files from multiple different individuals at the same time, we can replace the BAM file with  **-b bamlist.txt**, which would be a file we make in vi that is just a list of BAM files For example:

```
$ cat bamlist.txt
mappedReads_sample1.sorted.filtered.bam
mappedReads_sample2.sorted.filtered.bam
mappedReads_sample3.sorted.filtered.bam
mappedReads_sample4.sorted.filtered.bam
```

2.  `call` identifies positions where the sample differs from the reference genome.

![](https://github.com/nomascus/ANT3475_2026/blob/main/FILES/mpileupCall2.png?raw=true)

Here, you can see that there are two potential SNPs, a possible A/G at position 10, and a possible A/C at position 18.

With **bcftools call** we are taking the output of **bcftools mpileup** from the pipe, pulling out the variable positions, and converting the format to VCF. 

Here is what the other options from bcftools call mean: 

**-mv** use the standard variant caller and only take the variable positions
**-Ov** output format is VCF
**-o** name of the output file 

The output is a **VCF file** that only lists the variable postions that differ from the reference genome in one or more individuals. Here is an example of one.
```
##fileformat=VCFv4.2
##source=demo_pipeline
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1 sample2
chr1	105	.	A	G	50	PASS	.	GT	0/1	0/1
chr1	220	.	C	T	12	LowQual	.	GT	0/0	0/1
chr1	350	.	G	A	99	PASS	.	GT	1/1 1/0
```
As we discussed before, a VCF contains the following lines 

Metadata header lines that start with ##: 
```
##fileformat=VCFv4.2
##source=demo_pipeline
```

A column header line that starts with #:
`#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT sample1 sample2` 

and lines for each SNP 
```
chr1	105	.	A	G	50	PASS	.	GT	0/1	0/1
chr1	220	.	C	T	12	LowQual	.	GT	0/0	0/1
chr1	350	.	G	A	99	PASS	.	GT	1/1 1/0
```

# 5. Filter the VCF

The VCF we made includes any potential SNPs, even those that the program is not confident about. For example, if a nucleotide reference position is A, and it is covered by 9 reads 7 of which are A and two of which are T, it's difficult to be confident that this is a real SNP. In this case the confidence in the SNP call might be low.

To be safe, we should remove these low confidence SNPs

```
bcftools filter -i 'QUAL>30' HumanSamples.vcf -o HumanSamples.filtered.vcf` 
```

Here we are using another bcftools program called **bcftools filter**  to only keep SNPs with a quality score of at least 30. 

The option **-i** stands for "information" and allows us to write a complex string of filtering parameters. For now we are just going to tell it about the quality scores, but remember to include the quotes.


----------

## Summary of File Types and programs


| File Type | Description | Program to use|
| ----------- | ----- |----- |
| FASTQ | Sequencing reads| bwa mem |
| SAM | Mapping file | samtools view, samtools sort |
| BAM | Compressed mapping file | samtools filter,  bcftools mpileup, bcftools call|
| BAI | BAM index file | samtools index |
| .amb, .ann, .bwt, .fai, .pac, and .sa| reference genome index files| samtools index, bwa index|
| VCF | Variant call file | bcftools filter |


