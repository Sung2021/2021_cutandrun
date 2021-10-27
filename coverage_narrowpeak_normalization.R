library(ggplot2)
library(dplyr)
library(DESeq2)
library(GenomicRanges)
setwd('~/Desktop/HMH/')
theme_set(theme_bw())
## read input merged bed from macs2 call
dir <- 'rds/cutandrun/coverage'
list.files(dir, full.names = T)

### when reading coverage bed
for(i in 1:length(list.files(dir, full.names = T))){
  bed <- read.table(list.files(dir, full.names = T)[i])
  colnames(bed) <- c('chr','start','end',
                            'peak','score','strand',
                            'signalvalue','pvalue','qvalue','peak.score',
                            'read','frg_len','ref_len','coverage')
  assign(substr(list.files(dir, full.names = T)[i], 24,39),bed)
}
ls()
ls(pattern = 'effector') 
ls(pattern = 'Tcm') 
ls(pattern = 'Tem') 
ls(pattern = 'naive') 
ls(pattern = 'union')

eff1 <- effector_Tle3_1.
eff2 <- effector_Tle3_2.
nai1 <- naive_Tle3_1.36.
nai2 <- naive_Tle3_2.36.
nai3 <- naive_Tle3_3.36.
tcm1 <- Tcm_Tle3_1.36.q3
tcm2 <- Tcm_Tle3_2.36.q3
tem1 <- Tem_Tle3_1.36.q3
tem2 <- Tem_Tle3_2.36.q3

bed.list <- list(nai1,nai2,nai3,eff1,eff2,tcm1,tcm2,tem1,tem2)
names(bed.list) <- c('nai1','nai2','nai3','eff1','eff2','tcm1','tcm2','tem1','tem2')

for(i in 1:length(names(bed.list))){
  paste0('rds/cutandrun/',names(bed.list)[i],'.rds')
  saveRDS(bed.list[[i]], paste0('rds/cutandrun/',names(bed.list)[i],'.rds'))
}

for(i in 1:length(names(bed.list))){
  colnames(bed.list[[i]])[11] <- paste0(colnames(bed.list[[i]])[11],'_', names(bed.list)[i])
}

bed.list[[i]] %>% head()

signalvalue <- bed.list[[i]][,c(1:4)]
for(i in 1:9){
  signalvalue[,colnames(bed.list[[i]])[11]] <- bed.list[[i]][,11]
}

signalvalue %>% head()
saveRDS(signalvalue, 'rds/cutandrun/readvalue.2021.10.27.rds')

raw.value <- signalvalue[,c(5:13)]
raw.value %>% head()

group <- 'eff'
tmp <- raw.value[,grep(group, colnames(raw.value))]
tmp.avg <- tmp
tmp.avg[,'avg'] <- rowMeans(tmp.avg)
tmp.avg[,c(1:length(grep(group, colnames(tmp.avg))))]/mean(tmp.avg$avg)
assign(paste0(group, '_normalized'), tmp.avg[,c(1:length(grep(group, colnames(tmp.avg))))]/mean(tmp.avg$avg))

samples <- list()
for(group in c('nai','eff','tem','tcm')){
  tmp <- raw.value[,grep(group, colnames(raw.value))]
  tmp.avg <- tmp
  tmp.avg[,'avg'] <- rowMeans(tmp.avg)
  tmp.avg[,c(1:length(grep(group, colnames(tmp.avg))))]/mean(tmp.avg$avg)
  assign(paste0(group, '_normalized'), tmp.avg[,c(1:length(grep(group, colnames(tmp.avg))))]/mean(tmp.avg$avg))
  samples[['group']] <- tmp.avg[,c(1:length(grep(group, colnames(tmp.avg))))]/mean(tmp.avg$avg)
}
ls(pattern = 'normalized')

tmp.normalized <- signalvalue[,c(1:4)]
for(group in c('nai','eff','tem','tcm')){
  tmp <- raw.value[,grep(group, colnames(raw.value))]
  tmp.avg <- tmp
  tmp.avg[,'avg'] <- rowMeans(tmp.avg)
  tmp.avg[,c(1:length(grep(group, colnames(tmp.avg))))]/mean(tmp.avg$avg)
  tmp.normalized <- cbind(tmp.normalized, tmp.avg[,c(1:length(grep(group, colnames(tmp.avg))))]/mean(tmp.avg$avg))
}

tmp.normalized %>% head()
tmp.normalized %>% dim()
write.csv(tmp.normalized, 'rds/cutandrun/cutandrun.peak.normalized.csv')


normalization_factor <- list()
for(group in c('nai','eff','tem','tcm')){
  tmp <- raw.value[,grep(group, colnames(raw.value))]
  tmp.avg <- tmp
  tmp.avg[,'avg'] <- rowMeans(tmp.avg)
  normalization_factor[[group]] <- mean(tmp.avg$avg)
}
normalization_factor

