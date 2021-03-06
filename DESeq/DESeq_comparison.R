library(dplyr)
library(ggplot2)
library(DESeq2)
library(edgeR)
library(reshape)

####################################################################
########### DESeq2 analysis ###############
## read.value : raw count matrix
## check the duplication 
duplicated(read.value$peak) %>% table()
rownames(read.value) <- read.value$peak
## file location:  'rds/cutandrun/readvalue.2021.10.28.rename.rds')

## compare two groups : Tcm and Tem
count.mtx <- read.value[,grep('tcm|tem', colnames(read.value))]
group1 <- 'tcm' ## control
group2 <- 'tem' ## test

### generate the information sheet
info <- data.frame(matrix(nrow = length(colnames(count.mtx)), ncol = 3))
colnames(info) <- c('sample', 'cell_type','mutation')
info$sample <- substr(colnames(count.mtx),6,9)
info$cell_type <- substr(colnames(count.mtx),6,8)
info$mutation <- 'no'
info$cell_type <- factor(info$cell_type, levels = c(group1,group2))
info$mutation <- factor(info$mutation)
info

### generate the DESeq object
dds <- DESeqDataSetFromMatrix(count.mtx, info, ~cell_type)
dds <- DESeq(dds)

### vst for PCA : this is a quick PCA, for the normalized data_PCA, use prcomp (in other code) ### 
vsd <- vst(dds,blind=TRUE)
plotPCA(vsd, intgroup="cell_type") # to remove the dependence of the variance on the mean
pcaplot <- plotPCA(vsd, intgroup="sample", return=T)  #using the DESEQ2 plotPCA fxn we can

### drawing plot
theme <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), 
               axis.line = element_line(colour = "black"))
ggplot(pcaplot, aes(PC1,PC2, shape=sample)) + geom_point(size=3, alpha=0.8) + theme

### log2FC quick look ###
res <- results(dds)
set1 <- paste0('rds/cutandrun/res.', group1,'_',group2)
set2 <- paste0('rds/cutandrun/dds.', group1,'_',group2)
write.csv(res, paste0(set1, '.csv'))
saveRDS(dds, paste0(set2, '.rds'))

hist(res$log2FoldChange, breaks = 100)
hist(res$padj, breaks = 100)
## quick look for the volcano plot
EnhancedVolcano::EnhancedVolcano(res,
                                 lab = rownames(res),
                                 pCutoff = 0.05, FCcutoff = 1,
                                 x = 'log2FoldChange',
                                 y = 'padj',
                                 pointSize = 1.0,
                                 labSize = 4.0)

## calculate the number of peaks in each condition
res <- results(dds)
res <- data.frame(res)
res %>% filter(padj <= 0.05) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 2 | log2FoldChange <= -2) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1.5 | log2FoldChange <= -1.5) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1 | log2FoldChange <= -1) %>% dim()
