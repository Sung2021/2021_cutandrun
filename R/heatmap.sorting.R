### sorting the peaks for the heatmap visualization ###
norm.log2 <- read.csv('rds/cutandrun/11.09/cutandrun.2021.11.10.norm.log2FC.dif_peak_all.csv', row.names = 1,
                      stringsAsFactors = F)
norm.log2 %>% head()

table(norm.log2$dif_peak != 'NA' )
norm.log2.dif.peak <- norm.log2 %>% filter(dif_peak !='NA')
norm.log2.dif.peak %>% dim()
norm.log2.dif.peak %>% head()
write.csv(norm.log2.dif.peak, 'rds/cutandrun/11.09/cutandrun.2021.11.11.norm.log2FC.dif_peak_all.filtered.csv')

tmp <- norm.log2.dif.peak %>% select(contains('read_') | dif_peak) %>% arrange(dif_peak)
tmp %>% head()

## 2021.11.11 ##
tmp2 <- tmp
tmp2 <- tmp2[,-c(1:2)]
tmp2 %>% head()
dif_peaks <- factor(tmp2$dif_peak)
levels(dif_peaks)

tmp3 <- data.frame()
for(i in 1:14){
  input <- tmp2 %>% filter(dif_peaks == levels(dif_peaks)[i])
  input %>% head()
  input.avg <- data.frame(row.names = rownames(input))
  input.avg[,'nai'] <- input[,grep('nai', colnames(input))] %>% rowMeans()
  input.avg[,'eff'] <- input[,grep('eff', colnames(input))] %>% rowMeans()
  input.avg[,'tem'] <- input[,grep('tem', colnames(input))] %>% rowMeans()
  input.avg[,'tcm'] <- input[,grep('tcm', colnames(input))] %>% rowMeans()
  input.avg %>% head()
  row.order <- input.avg %>% arrange(eff,tem,tcm) %>% rownames()
  #row.order <- input.avg %>% arrange(desc(eff),desc(tem),desc(tcm)) %>% rownames()
  input <- input[row.order,]
  tmp3 <- rbind(tmp3,input)
}
df.anno <- data.frame(row.names = rownames(tmp3))
df.anno$group <- tmp3$dif_peak
x <- pheatmap::pheatmap(zscore(log2(tmp3[,-10]+1)), show_rownames = F, cluster_rows = F, cluster_cols = F,
                   col = colorRampPalette(c("navy", "white", "red"))(1000),
                   fontsize_row = 6,fontsize_col = 10,gaps_row = vec2,annotation_row = df.anno)




grep('A|B|C',tmp$dif_peak) %>% length()
grep('D|E|F',tmp$dif_peak) %>% length()
tmp$dif_peak %>% length()


pheatmap::pheatmap(tmp[,-12], show_rownames = F, cluster_rows = F, cluster_cols = F)
pheatmap::pheatmap(log2(tmp[,-12]+1), show_rownames = F, cluster_rows = F, cluster_cols = F)
## zscore function 
zscore <- function(input.data = input.data){
  input.data.rowsums <- rowSums(input.data)
  input.data.mean <- rowMeans(input.data)
  input.data.sd <- matrixStats::rowSds(as.matrix(input.data))
  names(input.data.sd) <- rownames(input.data)
  zscore <- (input.data-input.data.mean)/input.data.sd
  return(zscore)
}

pheatmap::pheatmap(zscore(log2(tmp[,-12]+1)), show_rownames = F, cluster_rows = F, cluster_cols = F,
                   col = colorRampPalette(c("navy", "white", "red"))(1000),
                   fontsize_row = 6,fontsize_col = 10)

df.anno <- data.frame(row.names = rownames(tmp))
df.anno$group <- tmp$dif_peak
table(df.anno$group)

pheatmap::pheatmap(zscore(log2(tmp[,-12]+1)), show_rownames = F, cluster_rows = F, cluster_cols = F,
                   col = colorRampPalette(c("navy", "white", "red"))(1000),
                   fontsize_row = 6,fontsize_col = 10, annotation_row = df.anno)


pheatmap::pheatmap(zscore(log2(tmp[,-c(12,1:2)]+1)), show_rownames = F, cluster_rows = F, cluster_cols = F,
                   col = colorRampPalette(c("navy", "white", "red"))(1000),
                   fontsize_row = 6,fontsize_col = 10, annotation_row = df.anno,
                   gaps_row = vec2)

vec <- c(632,530,689,70,180,57,84,
         202,442,435,8,234,45,19)
vec2 <- vector()

tmp.vec <- vector()
for(j in 1:length(vec)){
  print(j)
  print(vec[j-1])
  sum.tmp = sum(tmp.vec,vec[j-1])
  vec2[j] <- sum.tmp
  tmp.vec <- sum.tmp
}
vec2 <- vec2[-1] %>% as.vector()

x <- pheatmap::pheatmap(zscore(log2(tmp[,-c(12,1:2)]+1)), show_rownames = F, cluster_rows = F, cluster_cols = F,
                        col = colorRampPalette(c("navy", "white", "red"))(1000),
                        fontsize_row = 6,fontsize_col = 10, annotation_row = df.anno,
                        gaps_row = vec2)


save_pheatmap_pdf <- function(x, filename, width=7, height=12) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf(x, "heatmap.11.11.pdf", height = 30)


