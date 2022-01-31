#!/bin/bash

path=~/Desktop/HMH/rds/cutandrun/narrowpeaks
cd $path
for input in *Peak;
 do 
 bedtools coverage -a $path/union/union_peaks.narrowPeak -b ${input} > $path/coverage_bed/${input}.cov.bed
 done
 
