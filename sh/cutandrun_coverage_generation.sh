#!/bin/bash

### bam file to paired bed file for bedtools coverage input
path=~/Desktop/Sung_work/raw_tmp/TLE3_cut_run
cd $path/sorted_bam
#for input in *bam; do bedtools bamtobed -i $path/sorted_bam/${input} > $path/bed/${input}.pe.bed; done
cd $path/bed
for input in *bed; do bedtools coverage -a $path/macs2_all/union_peaks.narrowPeak -b $path/bed/${input} > $path/coverage/${input}.cov.bed; done
