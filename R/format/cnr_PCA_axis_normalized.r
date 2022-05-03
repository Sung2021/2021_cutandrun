## by Shaoqi
norm.peak.all <- read.csv('cutandrun/re_normalized.11.15/cutandrun.2021.11.17.cleared.peak.log2fc.csv', 
                          header = T, row.names = 1)
norm.peak.all %>% head()
norm.peak.all %>% dim() ## 34271
norm.peak <- norm.peak.all[,grep('read_', colnames(norm.peak.all))]
norm.peak %>% dim()
norm.peak <- norm.peak[,c(3:ncol(norm.peak))]
norm.peak %>% colnames()
norm.peak[1:3,]
pca = as.data.frame(prcomp(t(norm.peak), scale. = T)$x)
pca %>% str()

ggplot(pca, aes(PC1,PC2,label=rownames(pca), 
                shape=substr(rownames(pca), 6,8))) + geom_point()

PCA <- prcomp(t(norm.peak),scale.=T)
PCA <- (PCA$sdev/sum(PCA$sdev))[1:2]
PCA_percentage <- paste(round(PCA,4)*10^2,'%',sep = '')
