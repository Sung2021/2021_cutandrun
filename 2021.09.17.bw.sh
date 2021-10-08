#!/bin/bash

##### this step can be done after macs2 peak calling 
##### this step use pileup files from macs2 peak calling
##### only for visualization in IGV 
##### replicates are handled in for loop for the convenience

path=~/Desktop/kexin
## set input name of the replicates

### macs2 callpeak for pileup files (using no rmd files)
for input in con1 con2 eko1 eko2; do macs2 callpeak -t $path/bam/${input}.36.q30.sorted2.rmd.bam -f BAMPE -g mm -n ${input} -B -q 0.05 --SPMR --outdir $path/macs2 ; done

## sort the pileup files (from macs2 output) by chromosomal order
for input in con1 con2 eko1 eko2; do sort -k1,1 -k2,2n $path/macs2/${input}_treat_pileup.bdg > $path/macs2/${input}.sorted.bdg ; done

## covert the sorted pileup files to bw for the IGV use
for input in con1 con2 eko1 eko2; do bedGraphToBigWig $path/macs2/${input}.sorted.bdg ~/Desktop/software/mm10.chrom.sizes.txt $path/bw/${input}.sorted.bdg.bw ; done

## add delete intermediate files command below
aq