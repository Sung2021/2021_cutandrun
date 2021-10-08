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


input_folder=Tem_Tle3
input=Tem_Tle3_2
crop=36
path=~/Desktop/Sung_work/raw_tmp/TLE3_cut_run/${input_folder}
cd ~/Desktop/Sung_work/raw_tmp/TLE3_cut_run


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
-S $path/${input}_sam/${input}.${crop}.sam 2> $path/${input}_sam/${input}.${crop}.bowtie2.log


#############
### sam to bam conversion
### filter out MAPQ <30 
### sorting bam
### indexing bam
### remove intermediate files (sam)
### path again for the location of the sam

samtools view -@ 10 -S -b -@ 10 $path/${input}_sam/${input}.${crop}.sam > $path/${input}_sam/${input}.${crop}.bam
samtools view -@ 10 -q 30 -b $path/${input}_sam/${input}.${crop}.bam > $path/${input}_sam/${input}.${crop}.q30.bam

samtools sort -@ 10 $path/${input}_sam/${input}.${crop}.q30.bam -o $path/${input}_sam/${input}.${crop}.q30.sorted.bam
samtools index -@ 10 $path/${input}_sam/${input}.${crop}.q30.sorted.bam 

rm $path/${input}_sam/${input}.${crop}.sam

### remove duplicates 
### picard tool 
### output: bam with removal of duplicates

### picard markduplicates
java -jar ~/Desktop/software/picard.jar MarkDuplicates \
I= $path/${input}_sam/${input}.${crop}.q30.sorted.bam  \
O= $path/${input}_sam/${input}.${crop}.q30.sorted.rmd.bam \
M= $path/${input}_sam/${input}.${crop}.q30.sorted.rmd.metrics.txt \
ASSUME_SORTED=TRUE  REMOVE_DUPLICATES=true CREATE_INDEX=true QUIET=true 

### sort rmd bams
### output : sorted and duplicates-removed bam
samtools sort -@ 10 $path/${input}_sam/${input}.${crop}.q30.sorted.rmd.bam -o $path/${input}_sam/${input}.${crop}.q30.sorted.sorted.rmd.bam

################################################################
####### extra steps to count the reads through the flow ########

echo 'bam reads'
samtools view -@ 10 -c $path/${input}_sam/${input}.${crop}.bam
echo 'q30 sorted rmd bam sorted reads'
samtools view -@ 10 -c $path/${input}_sam/${input}.${crop}.q30.sorted.sorted.rmd.bam


