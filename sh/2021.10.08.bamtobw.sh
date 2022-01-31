#!/bin/bash


path=~/Desktop/Sung_work/raw_tmp/TLE3_cut_run
cd ~/Desktop/Sung_work/raw_tmp/TLE3_cut_run

input=naive_Tle3
### merge bams
samtools merge -@ 10 -f -o $path/${input}.merged.bam \
$path/${input}/${input}_1_sam/${input}_1.36.q30.sorted.sorted.rmd.bam $path/${input}/${input}_2_sam/${input}_2.36.q30.sorted.sorted.rmd.bam

### macs2 callpeak
macs2 callpeak -t $path/${input}.merged.bam -f BAMPE -g mm -n ${input} -B -q 0.05 --SPMR --outdir $path/macs2

### sort the pileup files
sort -k1,1 -k2,2n $path/macs2/${input}_treat_pileup.bdg > $path/macs2/${input}.sorted.bdg

### bedGraphToBigWig
bedGraphToBigWig $path/macs2/${input}.sorted.bdg ~/Desktop/software/mm10.chrom.sizes.txt $path/bw/${input}.sorted.bdg.bw
