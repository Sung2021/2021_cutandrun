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

tmp.normalized <- signalvalue[,c(1:4)]
for(group in c('nai','eff','tem','tcm')){
  tmp <- raw.value[,grep(group, colnames(raw.value))]
  tmp.avg <- tmp
  tmp.avg[,'avg'] <- rowMeans(tmp.avg)
  tmp.avg[,c(1:length(grep(group, colnames(tmp.avg))))]/mean(tmp.avg$avg)
  tmp.normalized <- cbind(tmp.normalized, 
                          tmp.avg[,c(1:length(grep(group, colnames(tmp.avg))))]/mean(tmp.avg$avg))
}
ls(pattern = 'normalized')
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


#########################################################################
############ 2021.10.28 
read.value <- readRDS('rds/cutandrun/readvalue.2021.10.27.rds')
read.value %>% head()

tmp.normalized <- read.value[,c(1:4)]
for(group in c('nai','eff','tem','tcm')){
  tmp <- read.value[,grep(group, colnames(read.value))]
  tmp.avg <- tmp
  tmp.avg[,'avg'] <- rowMeans(tmp.avg)
  tmp.avg[,c(1:length(grep(group, colnames(tmp.avg))))]/mean(tmp.avg$avg)
  tmp.normalized <- cbind(tmp.normalized, 
                          tmp.avg[,c(1:length(grep(group, colnames(tmp.avg))))]/mean(tmp.avg$avg))
}
ls(pattern = 'normalized')
tmp.normalized %>% head()
tmp.normalized %>% dim()

normalization_factor <- list()
for(group in c('nai','eff','tem','tcm')){
  tmp <- read.value[,grep(group, colnames(read.value))]
  tmp.avg <- tmp
  tmp.avg[,'avg'] <- rowMeans(tmp.avg)
  normalization_factor[[group]] <- mean(tmp.avg$avg)
}
normalization_factor
## nai 108.8927, eff 83.19871, tem 125.3606 tcm 96.11367

write.csv(read.value, 'rds/cutandrun/cutandrun.peak.raw.csv')



####################################################################
########### DESeq2 analysis ###############

read.value %>% head()
duplicated(read.value$peak) %>% table()
rownames(read.value) <- read.value$peak
saveRDS(read.value, 'rds/cutandrun/readvalue.2021.10.28.rename.rds')
read.value %>% colnames()



count.mtx <- read.value[,grep('tcm|tem', colnames(read.value))]
count.mtx %>% head()
group1 <- 'tcm'
group2 <- 'tem'


### information sheet
info <- data.frame(matrix(nrow = length(colnames(count.mtx)), ncol = 3))
colnames(info) <- c('sample', 'cell_type','mutation')
info$sample <- substr(colnames(count.mtx),6,9)
info$cell_type <- substr(colnames(count.mtx),6,8)
info$mutation <- 'no'
info
info$cell_type <- factor(info$cell_type, levels = c(group1,group2))
info$mutation <- factor(info$mutation)
str(info)
info

dds <- DESeqDataSetFromMatrix(count.mtx, info, ~cell_type)
dds <- DESeq(dds)

### vst for PCA ### 
vsd <- vst(dds,blind=TRUE)
plotPCA(vsd, intgroup="cell_type") # to remove the dependence of the variance on the mean
pcaplot <- plotPCA(vsd, intgroup="sample", return=T)  #using the DESEQ2 plotPCA fxn we can
### drawing plot
theme <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), 
               axis.line = element_line(colour = "black"))
ggplot(pcaplot, aes(PC1,PC2, shape=sample)) + geom_point(size=3, alpha=0.8) + 
  theme

### log2FC quick look ###
res <- results(dds)
set1 <- paste0('rds/cutandrun/res.', group1,'_',group2)
set2 <- paste0('rds/cutandrun/dds.', group1,'_',group2)
write.csv(res, paste0(set1, '.csv'))
saveRDS(dds, paste0(set2, '.rds'))

hist(res$log2FoldChange, breaks = 100)
hist(res$padj, breaks = 100)
EnhancedVolcano::EnhancedVolcano(res,
                                 lab = rownames(res),
                                 pCutoff = 0.05, FCcutoff = 1,
                                 x = 'log2FoldChange',
                                 y = 'padj',
                                 pointSize = 1.0,
                                 labSize = 4.0)

res <- results(dds)
res %>% dim()
res <- data.frame(res)
res %>% filter(padj <= 0.05) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 2 | log2FoldChange <= -2) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1.5 | log2FoldChange <= -1.5) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1 | log2FoldChange <= -1) %>% dim()



for(i in seq_along(list.files(dir, pattern = 'csv'))){
  dir <- 'rds/cutandrun/res'
  tmp <- read.csv(list.files(dir, pattern = 'csv', full.names = T)[i], header = T)
  colnames(tmp)[3] <- paste0(colnames(tmp)[3], '_',substr(list.files(dir, pattern = 'csv')[i],5,11) )
  colnames(tmp)[6] <- paste0(colnames(tmp)[6], '_',substr(list.files(dir, pattern = 'csv')[i],5,11) )
  colnames(tmp)[7] <- paste0(colnames(tmp)[7], '_',substr(list.files(dir, pattern = 'csv')[i],5,11) )
  assign(substr(list.files(dir, pattern = 'csv')[i],1,11), tmp)
}
ls(pattern = 'res.')
res.eff_tcm %>% head()

res.all <- data.frame(matrix(nrow=nrow(res.eff_tcm), ncol = 1))
rownames(res.all) <- paste0('union_peak_', rownames(res.all))
res.all %>% head()
for(i in seq_along(list.files(dir, pattern = 'csv'))){
  dir <- 'rds/cutandrun/res'
  tmp <- read.csv(list.files(dir, pattern = 'csv', full.names = T)[i], header = T)
  colnames(tmp)[3] <- paste0(colnames(tmp)[3], '_',substr(list.files(dir, pattern = 'csv')[i],5,11) )
  colnames(tmp)[6] <- paste0(colnames(tmp)[6], '_',substr(list.files(dir, pattern = 'csv')[i],5,11) )
  colnames(tmp)[7] <- paste0(colnames(tmp)[7], '_',substr(list.files(dir, pattern = 'csv')[i],5,11) )
  res.all <- cbind(res.all, tmp[,c(3,6,7)])
}
res.all <- res.all[,-1]
res.all %>% head()
write.csv(res.all, 'rds/cutandrun/res/res.all.csv')





norm.peak <- read.csv('rds/cutandrun/normalized/cutandrun.peak.normalized.csv', header = T)
rownames(norm.peak) <- norm.peak$peak

colnames(res.all)
#group <- 'nai_eff'
#group <- 'nai_tcm'
#group <- 'nai_tem'
#group <- 'eff_tcm'
#group <- 'eff_tem'
#group <- 'tcm_tem'

res.tmp <- res.all[,grep(group, colnames(res.all))]
for(i in 1:3){
  colnames(res.tmp)[i] <- substr(colnames(res.tmp)[i], 1, nchar(colnames(res.tmp)[i])-8)
}
res.tmp %>% head()

for(i in c(2,1.5,1)){
  threshold <- i
  res.tmp2 <- res.tmp %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= threshold | log2FoldChange <= -threshold)
  res.tmp2 <- cbind(res.tmp2, norm.peak[rownames(res.tmp2),c(6:14,2:4)])
  res.tmp2 %>% head()
  group.save <- paste0('_log2FC.',threshold)
  write.csv(res.tmp2, paste0('rds/cutandrun/res/',group, group.save, '.csv'))
}

####################################################################
########### 2021.11.01 ################


library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

readvalue <- readRDS('rds/cutandrun/count_raw/readvalue.2021.10.28.rename.rds')
readvalue %>% head()
df <- readvalue[,c(1:3)]
gr = makeGRangesFromDataFrame(df, keep.extra.columns=T)
peakAnno = ChIPseeker::annotatePeak(gr, tssRegion=c(-1000, 1000), 
                                    TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                                    annoDb="org.Mm.eg.db")
gr %>% head()
peakAnno %>% str()
peakAnno@anno %>% head()
peakAnno@anno %>% data.frame() %>% dim()
peakAnno@anno %>% data.frame() %>% head()
peakAnno@anno %>% names() %>% length()
peakAnno.df <- peakAnno@anno %>% data.frame(check.names = F)
rownames(peakAnno.df) <- names(peakAnno@anno)
peakAnno.df %>% head()
saveRDS(peakAnno.df, 'rds/cutandrun/dds/peak_annotation.tss1000.total.2021.11.01.rds')




peakAnno = ChIPseeker::annotatePeak(gr, tssRegion=c(-3000, 3000), 
                                    TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                                    annoDb="org.Mm.eg.db")
gr %>% head()
peakAnno %>% str()
peakAnno@anno %>% head()
peakAnno@anno %>% data.frame() %>% dim()
peakAnno@anno %>% data.frame() %>% head()
peakAnno@anno %>% names() %>% length()
peakAnno.df <- peakAnno@anno %>% data.frame(check.names = F)
rownames(peakAnno.df) <- names(peakAnno@anno)

write.csv(peakAnno.df, 'rds/cutandrun/cutandrun.peak.annotation.csv')
