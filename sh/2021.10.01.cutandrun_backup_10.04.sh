#!/bin/bash

###### these codes will generate one final bam file from two fastqs (paired) ###
###### assign one input name for the all downstream process
###### in this case, trimming process remaining 36bp was used so 'crop' indicate the remaining sequence number as an indicator
### two variables 
### input, crop
### set path

##### general workflow #####
### Trimming of fastqs remaining {crop}bp 
### There are other trimming methods
### Alignment of samples : bowtie2 for cut and run
### output log from bowtie2 would be used to calculate mappable reads from fastq files.
### calculate only aligned reads in output log
### sam to bam conversion : samtools
### filter out MAPQ <30 : samtools
### sorting bam : samtools
### indexing bam : samtools
### remove intermediate files (sam)
### remove duplicates by picard markduplicates
### sort rmd bams : samtools
### extra steps to count the reads through the flow

### input folders : naive_Tle3 naive_IgG effector_Tle3 Tem_Tle3 Tcm_Tle3 


input_folder=naive_Tle3
input=naive_Tle3_1
crop=36
path=~/Desktop/Sung_work/raw_tmp/TLE3_cut_run/${input_folder}
cd ~/Desktop/Sung_work/raw_tmp/TLE3_cut_run

### Trimming of fastqs remaining {crop}bp 
### Trimmomatic program location : ~/Desktop/software/Trimmomatic-0.39/trimmomatic-0.39.jar
### output : trimmed fastqs (R1 and R2)
## run
java -jar ~/Desktop/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 10 \
$path/${input}.R1.fastq.gz $path/${input}.R2.fastq.gz \
$path/${input}.${crop}.R1.paired.fastq.gz $path/${input}.${crop}.R1.unpaired.fastq.gz \
$path/${input}.${crop}.R2.paired.fastq.gz $path/${input}.${crop}.R2.unpaired.fastq.gz \
CROP:${crop}

### Alignment of samples : bowtie2 for ATAC-seq
### bowtie2
### bowtie2 index location : ~/Desktop/ref/mm10/mm10.index
### indicate location for the output log from bowtie2 : ~/Desktop/Sung_work/log/
### indicate the location for the output sam : ~/Desktop/Sung_work/sam/
### output: bam

mkdir $path/${input}_sam

bowtie2 --mm \
-p 10 \
--no-unal \
--non-deterministic \
--no-discordant \
-x ~/Desktop/ref/mm10/mm10.index \
-1 $path/${input}.${crop}.R1.paired.fastq.gz \
-2 $path/${input}.${crop}.R2.paired.fastq.gz \
-S $path/${input}_sam/${input}.${crop}.sam 2> $path/${input}_sam/${input}.${crop}.log
