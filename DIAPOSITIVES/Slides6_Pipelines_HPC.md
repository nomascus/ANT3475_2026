# ANT3475 Class 6:  Pipelines and HPC

Here are the files to download for today's class. Please start my making a new directory as follows.
```
mkdir ~/scratch/6_HPC
cd ~/scratch/6_HPC

wget https://raw.githubusercontent.com/nomascus/ANT3475_2026/refs/heads/main/FILES/6_HPC/reference_mt.fasta
wget https://raw.githubusercontent.com/nomascus/ANT3475_2026/refs/heads/main/FILES/6_HPC/mt_multiSNP_2000_R1.fastq 
wget https://raw.githubusercontent.com/nomascus/ANT3475_2026/refs/heads/main/FILES/6_HPC/mt_multiSNP_2000_R2.fastq
wget https://raw.githubusercontent.com/nomascus/ANT3475_2026/refs/heads/main/FILES/6_HPC/reads_sample1_F.fastq
wget https://raw.githubusercontent.com/nomascus/ANT3475_2026/refs/heads/main/FILES/6_HPC/reads_sample1_R.fastq
wget https://raw.githubusercontent.com/nomascus/ANT3475_2026/refs/heads/main/FILES/6_HPC/metadataHPC.txt
wget https://raw.githubusercontent.com/nomascus/ANT3475_2026/refs/heads/main/FILES/6_HPC/metadata_batch.txt
wget https://raw.githubusercontent.com/nomascus/ANT3475_2026/refs/heads/main/FILES/6_HPC/mtDNA_mappingVCF.sh
wget https://raw.githubusercontent.com/nomascus/ANT3475_2026/refs/heads/main/FILES/6_HPC/mtDNA_mappingVCF2.sh
wget https://raw.githubusercontent.com/nomascus/ANT3475_2026/refs/heads/main/FILES/6_HPC/mtDNA_mappingVCF2.sbatch
```

# Putting it all together into a pipeline

We did a whole bunch of bioinformatics last week when we made those VCFs! However, it's not very efficient to type all these lines of text one by one. Especially if you want to run this on tens or hundreds samples. This is why we make scripts called **pipelines** that tie together all these programs one after another so that we can run the whole process in a single line of code.

A pipeline will be a long text file we make in vi that ends with .sh if we are running it on the command line  or .sbatch if we are running it remotely (more on that later).

Here is a an example pipeline for everything we did last week.

``` vi mtDNA_mappingVCF.sh ```

```
#!/bin/bash

############################
# 0. Load modules
############################

module load bwa
module load samtools
module load bcftools

############################
# 1. Define input files.
############################

REF=reference_mt.fasta
R1=reads_sample1_F.fastq
R2=reads_sample1_R.fastq
SAMPLE='Sample1'

############################
# 2. Index reference files 
############################

echo "Indexing reference..."

samtools faidx $REF
bwa index $REF

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
```



1) The first line is always **#!/bin/bash** This is a special line of code that test the computer this is a bash program. Don't forget it!

2) For all other lines a **#** means that the line is a comment line. Comment lines are **not** run by the computer (there are a few exceptions like #! and #SBATCH). You can use them make notes, make your code more readable and break it into visible blocks.

3) In block 0, we load the three programs, samtools, bwa, and bcftools 

4) In block 1,  We are telling the program the input files we want to use. Here we are declaring four variables, one for the reference assembly, one for the forward reads, one for the reverse reads, and one for the sample name. This is helpful, because if we want to run the pipeline on a different set of samples, references, or sets of reads we only need to change the names of these files one time.

5) In the remaining blocks, were are just calling the same commands as last week to map the reads and process the VCF, except we are mixing in the variable names we declared above. We've also added a few `echo` lines to print information to the screen so we can see what step is running.

To run the entire pipeline that maps reads and calls VCFs for this set of sequencing reads, all we need to do is type 

```bash mtDNA_mappingVCF.sh ```

and the entire pipeline runs!

### HPC Part 2

Remember when we discussed the environmental variables like $1, $2, $3 ? I told you not to use these, because they have a special function. If we set a variable to $1 in the pipeline, it takes the first item typed after the program and stores it as $1. Likewise, if we set a variable to $2 in the pipeline, it takes the second item typed after the program and stores it as $2. and so on and so forth.

What happens if we change block 1 of our pipeline to the following:
``` vi mtDNA_mappingVCF2.sh```
```
############################
# 1. Define input files.
############################

REF=$1 
R1=$2
R2=$3
SAMPLE=$4
```
When we call the program ```bash mtDNA_mappingVCF2.sh``` it  brings up a bunch of errors, because it's expecting four more pieces of information (in this case the reference, two sets of fastq reads, and a sample name), which it needs to run in the program

So if we run 
``` 
bash mtDNA_mappingVCF2.sh reference_mt.fasta reads_sample1_F.fastq reads_sample1_R.fastq Sample1
```
bash will store those four strings of text as $1, $2, $3, and $4 for use in the our script. To make things more comprehensible, we can then store those variables in more clear variables as follows

```
REF=$1 
R1=$2
R2=$3
SAMPLE=$4
```

| Environmental variable |Stored information |Variable name in pipeline
| - | - | - | - |
| $1 |reference_mt.fasta | $REF|
| $2 |reads_sample1_F.fastq | $R1|
| $3 |reads_sample1_R.fastq | $R2|
| $4 |Sample1 | $SAMPLE|

This grants the pipeline lots of flexibility, because we can just change the information on the command line and run the same pipeline! This becomes super useful if we want to run the pipeline on multiple samples or sets of reads at the same time.
We can do this with a **while** loop and a metadata file
```
[orkin@login1 6_Unix]$ cat metadataHPC.txt
reference_mt.fasta  reads_sample1_F.fastq  reads_sample1_R.fastq  Sample1
reference_mt.fasta  mt_multiSNP_2000_R1.fastq  mt_multiSNP_2000_R2.fastq  Sample1
```
```
while read line; do bash mtDNA_mappingVCF2.sh $line; done <metadataHPC.txt
```

Et voila! With this one line of coded we've run the gone from fastq to filtered VCF for two samples!

## What is an HPC

Until now, all of our work has taken place on the cloud computer at the command line. While the command line interface we've been learning is different than the graphical user interphase you are used to on Windows or Mac, you are probably wondering why all of this is necessary. Well, one of the main reasons we use the command line interface is so that we can run many jobs at the same time. Unlike with a regular computer, we can schedule tens, hundreds, or thousands of jobs to run simultaneously with only a few lines of code. 

A high-performance computing (HPC) cluster is a group of many computers connected together and managed so they can run large computational jobs efficiently. Instead of running analyses on your laptop, you submit jobs to the cluster, which runs them on powerful remote machines.

A typical cluster has three parts:

1) **The login node**. This is where we have been working. It's where you connect via ssh, edit scripts, move files, and submit jobs with `sbatch`. We don't use this for heavy computational tasks, because it's a shared space

2) **Compute nodes**. These are the powerful computers where we will send our scripts asking them to run our code. You don't ever log in to these nodes. 

3) **The job scheduler (SLURM)**. SLURM is a program that manages all the jobs being sent out to compute nodes by all the users. Think of it like a host at a restaurant with a long line of customers. It's the host's job to seat people at tables for the right amount of time and number of seats.

Think of this like a busy restaurant. It's the job of the host (SLURM) to seat a long line of customers (the jobs waiting in the queue) at empty tables of different sizes (the compute nodes) so they can eat.

To send our job out to the compute nodes, we need to add the following lines of code to the top of our script (just under #!/bin/bash)

```
#SBATCH --job-name=Joe_vcf_pipeline
#SBATCH --cpus-per-task=1  
#SBATCH --mem=1G  
#SBATCH --time=00:01:00  
#SBATCH --output=%x_%j.out  
#SBATCH --error=%x_%j.err
```
Let's break this down. The  **`#SBATCH`** lines tell the cluster what resources your job needs before it runs.

`#SBATCH --job-name=bwa_pipeline`
This gives your job a human-readable name, which makes output files easier to interpret. You can change this to whatever you want

`#SBATCH --cpus-per-task=1`
Requests the number of CPU cores for the job.  More CPUs result in a faster runtime for the job

`#SBATCH --mem=1G`
Requests the amount of RAM for the job. Too little ram and your job will crash, but too much ram and it will sit in the queue waiting to run

`#SBATCH --time=00:10:00`
Sets the maximum time requested for your job to run. Longer time requests sit in the queue for a while, but if your job is not finished before the requested time, it will die. Time is in hours:minutes:seconds

`#SBATCH --output=%x_%j.out`
This specifies a file name for the standard output (what would normally be printed to the screen). %x is code for the job name, and %j is for the job ID, a unique number the computer assigns to each job. For example, this might create a file called Joe_bwa_pipeline_12.out

`#SBATCH --error=%x_%j.err`
This specifies a file name for the standard error (error messages normally printed to the screen). This might create a file called Joe_bwa_pipeline_12.err


**In this pipeline:** it looks like this:

```vi  mtDNA_mappingVCF2.sbatch```
```
#!/bin/bash  
#SBATCH --job-name=Joe_vcf_pipeline  
#SBATCH --cpus-per-task=1  
#SBATCH --mem=1G  
#SBATCH --time=00:10:00  
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
```


To send the job to the compute nodes, we use the command **sbatch** instead of bash, and we can save the pipeline file with a different ending so we can identify it. Once a job is submitted, we are told by the cluster
```
[orkin@login1 TMP]$ sbatch mtDNA_mappingVCF2.sbatch
Submitted batch job 23
```
We can then check on the status of the job with the **squeue** command

```
[orkin@login1 TMP]$ squeue
JOBID USER  ACCOUNT NAME  ST  TIME_LEFT NODES CPUS TRES_PER_N MIN_MEM NODELIST (REASON)
23  orkin def-sponsor0 bwa_pipeline  R 0:58 1  1  N/A  1G nodecpupool1 (None)
```
This tells us info about each job in the queue. Most importantly, who is running is, the name of the job, the status, and the time left. Note the statuses **PD** = waiting in queue; **CF** = starting up, **R** = running

If the cluster has been idle for a while, you might see your job stuck in CF status for a few minutes. Once the first job runs the rest will run quickly. If your jobs get stuck and you want to remove them all from the queue, you can type **scancel** followed by your userid. For example: ```scancel orkin```

And note that we can run the pipeline on an enormous number samples at the same time using  a while loop with a metadata file!

Each line in the metadata file corresponds to a tab separated list of the necessary text $1 $2 $3 $4 to run the pipeline
```
[orkin@login1 TMP]$ cat metadata_batch.txt
reference_mt.fasta reads_sample1_F.fastq  reads_sample1_R.fastq  Sample1
reference_mt.fasta reads_sample1_F.fastq  reads_sample1_R.fastq  Sample2
reference_mt.fasta reads_sample1_F.fastq  reads_sample1_R.fastq  Sample3
reference_mt.fasta reads_sample1_F.fastq  reads_sample1_R.fastq  Sample4
reference_mt.fasta  mt_multiSNP_2000_R1.fastq  mt_multiSNP_2000_R2.fastq  Sample5
reference_mt.fasta  mt_multiSNP_2000_R1.fastq  mt_multiSNP_2000_R2.fastq  Sample6
reference_mt.fasta  mt_multiSNP_2000_R1.fastq  mt_multiSNP_2000_R2.fastq  Sample7
reference_mt.fasta  mt_multiSNP_2000_R1.fastq  mt_multiSNP_2000_R2.fastq  Sample8
```

```
while read line; do sbatch mtDNA_mappingVCF2.sbatch $line; done < metadata_batch.txt
```

And just like that, we've run the pipeline 8 times on different samples! 
