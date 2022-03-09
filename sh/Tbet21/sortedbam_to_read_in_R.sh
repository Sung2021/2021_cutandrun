## 1.merge effector Tle3 bam files
samtools merge Tle3.eff.merged.bam effector_Tle3_1.36.q30.sorted.sorted.rmd.bam effector_Tle3_2.36.q30.sorted.sorted.rmd.bam 
## 2.call peaks of merged Tle3 effector
macs2 callpeak --cutoff-analysis -t Tle3.eff.merged.bam -f BAMPE -g mm -n Tle3.eff -B -q 0.05 --SPMR --outdir macs2
## 3.bamtobed for individual Tle3 effector files
for input in *.sorted.rmd.bam;do bedtools bamtobed -i ${input} > ${input}.bamtobed.bed;done
## 4.calculate coverage using merged narrowPeak as a control 
for input in *.bamtobed.bed;do bedtools coverage -a macs2/Tle3.eff_peaks.narrowPeak -b ${input} > coverage/${input}.cov.bed;done
## 5.read files in R
tle3.r1 <- read.table('rds/cutandrun/Tle3.effector/coverage/effector_Tle3_1.36.q30.sorted.sorted.rmd.bam.bamtobed.bed.cov.bed')
tle3.r2 <- read.table('rds/cutandrun/Tle3.effector/coverage/effector_Tle3_2.36.q30.sorted.sorted.rmd.bam.bamtobed.bed.cov.bed')
colnames.bed <- c('chr','start','end',
                  'peak','score','strand',
                  'signalvalue','pvalue','qvalue','peak.score',
                  'read','frg_len','ref_len','coverage')
colnames(tle3.r1) <- colnames.bed
colnames(tle3.r2) <- colnames.bed
