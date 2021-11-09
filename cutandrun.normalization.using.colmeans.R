count.mtx <- read.value[,grep('nai|eff', colnames(read.value))]
group1 <- 'nai'
group2 <- 'eff'
count.mtx %>% head()


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

res <- results(dds)
res %>% dim()
res <- data.frame(res)
res %>% filter(padj <= 0.05) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 2 | log2FoldChange <= -2) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1.5 | log2FoldChange <= -1.5) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1 | log2FoldChange <= -1) %>% dim()

res %>% head()
res <- res[,c(2,5,6)]
rownames(res) <- rownames(norm.peak)
write.csv(res,paste0('rds/cutandrun/res.11.09/', group1,'_',group2, '.csv'))
res %>% head()
res.tmp <- res

for(i in c(2,1.5,1)){
  threshold <- i
  res.tmp2 <- res.tmp %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= threshold | log2FoldChange <= -threshold)
  res.tmp2 <- cbind(res.tmp2, norm.peak[rownames(res.tmp2),c(6:14,2:4)])
  colnames(res.tmp2) <- paste0(group1,'_', group2, '_',colnames(res.tmp2))
  res.tmp2 %>% head()
  group.save <- paste0('_log2FC.',threshold)
  write.csv(res.tmp2, paste0('rds/cutandrun/res.11.09/',paste0(group1,'_', group2, '_'), group.save, '.csv'))
}



##############################################################
##############################################################


count.mtx <- read.value[,grep('nai|tcm', colnames(read.value))]
group1 <- 'nai'
group2 <- 'tcm'
count.mtx %>% head()


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

res <- results(dds)
res %>% dim()
res <- data.frame(res)
res %>% filter(padj <= 0.05) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 2 | log2FoldChange <= -2) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1.5 | log2FoldChange <= -1.5) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1 | log2FoldChange <= -1) %>% dim()

res %>% head()
res <- res[,c(2,5,6)]
rownames(res) <- rownames(norm.peak)
write.csv(res,paste0('rds/cutandrun/res.11.09/', group1,'_',group2, '.csv'))
res %>% head()
res.tmp <- res

for(i in c(2,1.5,1)){
  threshold <- i
  res.tmp2 <- res.tmp %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= threshold | log2FoldChange <= -threshold)
  res.tmp2 <- cbind(res.tmp2, norm.peak[rownames(res.tmp2),c(6:14,2:4)])
  colnames(res.tmp2) <- paste0(group1,'_', group2, '_',colnames(res.tmp2))
  res.tmp2 %>% head()
  group.save <- paste0('_log2FC.',threshold)
  write.csv(res.tmp2, paste0('rds/cutandrun/res.11.09/',paste0(group1,'_', group2, '_'), group.save, '.csv'))
}



##############################################################
##############################################################


count.mtx <- read.value[,grep('nai|tem', colnames(read.value))]
group1 <- 'nai'
group2 <- 'tem'
count.mtx %>% head()


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

res <- results(dds)
res %>% dim()
res <- data.frame(res)
res %>% filter(padj <= 0.05) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 2 | log2FoldChange <= -2) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1.5 | log2FoldChange <= -1.5) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1 | log2FoldChange <= -1) %>% dim()

res %>% head()
res <- res[,c(2,5,6)]
rownames(res) <- rownames(norm.peak)
write.csv(res,paste0('rds/cutandrun/res.11.09/', group1,'_',group2, '.csv'))
res %>% head()
res.tmp <- res

for(i in c(2,1.5,1)){
  threshold <- i
  res.tmp2 <- res.tmp %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= threshold | log2FoldChange <= -threshold)
  res.tmp2 <- cbind(res.tmp2, norm.peak[rownames(res.tmp2),c(6:14,2:4)])
  colnames(res.tmp2) <- paste0(group1,'_', group2, '_',colnames(res.tmp2))
  res.tmp2 %>% head()
  group.save <- paste0('_log2FC.',threshold)
  write.csv(res.tmp2, paste0('rds/cutandrun/res.11.09/',paste0(group1,'_', group2, '_'), group.save, '.csv'))
}



##############################################################
##############################################################


count.mtx <- read.value[,grep('eff|tcm', colnames(read.value))]
group1 <- 'eff'
group2 <- 'tcm'
count.mtx %>% head()


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

res <- results(dds)
res %>% dim()
res <- data.frame(res)
res %>% filter(padj <= 0.05) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 2 | log2FoldChange <= -2) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1.5 | log2FoldChange <= -1.5) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1 | log2FoldChange <= -1) %>% dim()

res %>% head()
res <- res[,c(2,5,6)]
rownames(res) <- rownames(norm.peak)
write.csv(res,paste0('rds/cutandrun/res.11.09/', group1,'_',group2, '.csv'))
res %>% head()
res.tmp <- res

for(i in c(2,1.5,1)){
  threshold <- i
  res.tmp2 <- res.tmp %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= threshold | log2FoldChange <= -threshold)
  res.tmp2 <- cbind(res.tmp2, norm.peak[rownames(res.tmp2),c(6:14,2:4)])
  colnames(res.tmp2) <- paste0(group1,'_', group2, '_',colnames(res.tmp2))
  res.tmp2 %>% head()
  group.save <- paste0('_log2FC.',threshold)
  write.csv(res.tmp2, paste0('rds/cutandrun/res.11.09/',paste0(group1,'_', group2, '_'), group.save, '.csv'))
}



##############################################################
##############################################################


count.mtx <- read.value[,grep('eff|tem', colnames(read.value))]
group1 <- 'eff'
group2 <- 'tem'
count.mtx %>% head()


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

res <- results(dds)
res %>% dim()
res <- data.frame(res)
res %>% filter(padj <= 0.05) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 2 | log2FoldChange <= -2) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1.5 | log2FoldChange <= -1.5) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1 | log2FoldChange <= -1) %>% dim()

res %>% head()
res <- res[,c(2,5,6)]
rownames(res) <- rownames(norm.peak)
write.csv(res,paste0('rds/cutandrun/res.11.09/', group1,'_',group2, '.csv'))
res %>% head()
res.tmp <- res

for(i in c(2,1.5,1)){
  threshold <- i
  res.tmp2 <- res.tmp %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= threshold | log2FoldChange <= -threshold)
  res.tmp2 <- cbind(res.tmp2, norm.peak[rownames(res.tmp2),c(6:14,2:4)])
  colnames(res.tmp2) <- paste0(group1,'_', group2, '_',colnames(res.tmp2))
  res.tmp2 %>% head()
  group.save <- paste0('_log2FC.',threshold)
  write.csv(res.tmp2, paste0('rds/cutandrun/res.11.09/',paste0(group1,'_', group2, '_'), group.save, '.csv'))
}





##############################################################
##############################################################


count.mtx <- read.value[,grep('tcm|tem', colnames(read.value))]
group1 <- 'tcm'
group2 <- 'tem'
count.mtx %>% head()


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

res <- results(dds)
res %>% dim()
res <- data.frame(res)
res %>% filter(padj <= 0.05) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 2 | log2FoldChange <= -2) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1.5 | log2FoldChange <= -1.5) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1 | log2FoldChange <= -1) %>% dim()

res %>% head()
res <- res[,c(2,5,6)]
rownames(res) <- rownames(norm.peak)
write.csv(res,paste0('rds/cutandrun/res.11.09/', group1,'_',group2, '.csv'))
res %>% head()
res.tmp <- res

for(i in c(2,1.5,1)){
  threshold <- i
  res.tmp2 <- res.tmp %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= threshold | log2FoldChange <= -threshold)
  res.tmp2 <- cbind(res.tmp2, norm.peak[rownames(res.tmp2),c(6:14,2:4)])
  colnames(res.tmp2) <- paste0(group1,'_', group2, '_',colnames(res.tmp2))
  res.tmp2 %>% head()
  group.save <- paste0('_log2FC.',threshold)
  write.csv(res.tmp2, paste0('rds/cutandrun/res.11.09/',paste0(group1,'_', group2, '_'), group.save, '.csv'))
}





