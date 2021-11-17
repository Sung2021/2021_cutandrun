#!/bin/bash

path=~/Desktop/HMH/rds/cutandrun/bed/

mkdir $path/homer
cd $path

perl /Users/sung/opt/miniconda3/share/homer-4.10-0/bin/findMotifsGenome.pl 
for input in *bed; do perl /Users/sung/opt/miniconda3/share/homer-4.10-0/bin/findMotifsGenome.pl $path/${input} mm10 $path/homer/${input} ; done 
