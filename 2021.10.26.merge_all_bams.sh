#!/bin/bash


path=~/Desktop/Sung_work/raw_tmp/TLE3_cut_run
cd ~/Desktop/Sung_work/raw_tmp/TLE3_cut_run/sorted_bam

# samtools merge -@ 10 -f -o $path/all.merged.bam *bam 

# mkdir $path/macs2_all
# macs2 callpeak -t $path/all.merged.bam -f BAMPE -g mm -n union -B -q 0.05 --SPMR --outdir $path/macs2_all

#mkdir $path/macs2_merged_control
#for input in *bam; 
#do macs2 callpeak -t ${input} -f BAMPE -g mm -n ${input} -B -q 0.05 --SPMR --outdir $path/macs2_merged_control
#done 

