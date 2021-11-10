## cutandrun_2021.11.10.dif_peak.manipulation.and.save.log.R

norm.log2 <- read.csv('rds/cutandrun/11.09/cutandrun.2021.11.10.norm.log2FC.csv', row.names = 4)
norm.log2 %>% head()
norm.log2 %>% dim()
group <- 'nai_eff'
group2 <- 'read_nai'
#grep(group,colnames(norm.log2))
res <- norm.log2[,grep(group,colnames(norm.log2))]
res %>% head()
res %>% dim()
res %>% filter(.data[[colnames(res)[1]]] <= -log2(3)) %>% dim()
res <- res %>% filter(.data[[colnames(res)[1]]] <= -log2(3))
res <- res %>% filter(.data[[colnames(res)[3]]] < 0.1)
res %>% dim()
# norm.log2 %>% select(contains(group2)) %>% head()
res.norm <- norm.log2 %>% select(contains(group2))
res.norm <- res.norm[rownames(res),]
res.norm$group.avg <-  res.norm %>% rowMeans()
res.norm <- res.norm %>% filter(group.avg >= 0.7)
res.norm %>% dim()
table(res.norm$group.avg >= 0.7)
res.norm <- cbind(res.norm, norm.log2[rownames(res.norm),])
res.norm %>% head()
write.csv(res.norm, 'rds/cutandrun/Tn_Teff.dn.1087.csv') 
res.norm.eff.dn <- res.norm


group <- 'nai_tem'
group2 <- 'read_nai'
res <- norm.log2[,grep(group,colnames(norm.log2))]
res %>% head()
res %>% dim()
res %>% filter(.data[[colnames(res)[1]]] <= -log2(3)) %>% dim()
res <- res %>% filter(.data[[colnames(res)[1]]] <= -log2(3))
res <- res %>% filter(.data[[colnames(res)[3]]] < 0.1)
res %>% dim()
# norm.log2 %>% select(contains(group2)) %>% head()
res.norm <- norm.log2 %>% select(contains(group2))
res.norm <- res.norm[rownames(res),]
res.norm$group.avg <-  res.norm %>% rowMeans()
res.norm <- res.norm %>% filter(group.avg >= 0.7)
res.norm %>% dim()
table(res.norm$group.avg >= 0.7)
res.norm <- cbind(res.norm, norm.log2[rownames(res.norm),])
res.norm %>% head()
write.csv(res.norm, 'rds/cutandrun/Tn_Tem.dn.1156.csv') 
res.norm.tem.dn <- res.norm



group <- 'nai_tcm'
group2 <- 'read_nai'
res <- norm.log2[,grep(group,colnames(norm.log2))]
res %>% head()
res %>% dim()
res %>% filter(.data[[colnames(res)[1]]] <= -log2(3)) %>% dim()
res <- res %>% filter(.data[[colnames(res)[1]]] <= -log2(3))
res <- res %>% filter(.data[[colnames(res)[3]]] < 0.1)
res %>% dim()
# norm.log2 %>% select(contains(group2)) %>% head()
res.norm <- norm.log2 %>% select(contains(group2))
res.norm <- res.norm[rownames(res),]
res.norm$group.avg <-  res.norm %>% rowMeans()
res.norm <- res.norm %>% filter(group.avg >= 0.7)
res.norm %>% dim()
table(res.norm$group.avg >= 0.7)
res.norm <- cbind(res.norm, norm.log2[rownames(res.norm),])
res.norm %>% head()
write.csv(res.norm, 'rds/cutandrun/Tn_Tcm.dn.507.csv') 
res.norm.tcm.dn <- res.norm


ls(pattern = 'dn')

res.norm.eff.dn
res.norm.tem.dn
res.norm.tcm.dn

dn.eff.tem <- intersect(rownames(res.norm.eff.dn), rownames(res.norm.tem.dn))
dn.eff.tcm <- intersect(rownames(res.norm.eff.dn), rownames(res.norm.tcm.dn))
dn.tem.tcm <- intersect(rownames(res.norm.tem.dn), rownames(res.norm.tcm.dn))

dn.eff.tcm %>% length()

DEF_435 <- intersect(dn.eff.tem, dn.eff.tcm)
DE_442 <- setdiff(dn.eff.tem,DEF_435)
DF_8 <- setdiff(dn.eff.tcm,DEF_435)
EF_45 <- setdiff(dn.tem.tcm,DEF_435)
D_202 <- setdiff(setdiff(setdiff(rownames(res.norm.eff.dn), DE_442), DEF_435), DF_8)
E_234 <- setdiff(setdiff(setdiff(rownames(res.norm.tem.dn), DEF_435), DE_442), EF_45)
F_19 <- setdiff(setdiff(setdiff(rownames(res.norm.tcm.dn), DEF_435), EF_45), DF_8)


table(norm.log2$dif_peak)
norm.log2[DEF_435,'dif_peak'] <- 'DEF_435'
norm.log2[DE_442,'dif_peak'] <- 'DE_442'
norm.log2[DF_8,'dif_peak'] <- 'DF_8'
norm.log2[EF_45,'dif_peak'] <- 'EF_45'
norm.log2[D_202,'dif_peak'] <- 'D_202'
norm.log2[E_234,'dif_peak'] <- 'E_234'
norm.log2[F_19,'dif_peak'] <- 'F_19'

table(norm.log2$dif_peak) %>% data.frame()
norm.log2 %>% head()


write.csv(norm.log2, 'rds/cutandrun/11.09/cutandrun.2021.11.10.norm.log2FC.dif_peak_dn.csv')
write.csv(norm.log2, 'rds/cutandrun/11.09/cutandrun.2021.11.10.norm.log2FC.dif_peak_all.csv')

table(norm.log2$dif_peak != 'NA' )
norm.log2.dif.peak <- norm.log2 %>% filter(dif_peak !='NA')
norm.log2.dif.peak %>% dim()
norm.log2.dif.peak %>% head()
tmp <- norm.log2.dif.peak %>% select(contains('read_') | dif_peak) %>% arrange(dif_peak)
tmp %>% head()


grep('A|B|C',tmp$dif_peak) %>% length()
grep('D|E|F',tmp$dif_peak) %>% length()
tmp$dif_peak %>% length()
write.csv(norm.log2[rownames(tmp[grep('A|B|C',tmp$dif_peak), ]),],
          'rds/cutandrun/11.09/cutandrun.2021.11.10.norm.log2FC.dif_peak_up.csv' )
write.csv(norm.log2[rownames(tmp[grep('D|E|F',tmp$dif_peak), ]),],
          'rds/cutandrun/11.09/cutandrun.2021.11.10.norm.log2FC.dif_peak_dn.csv' )


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
                   col = colorRampPalette(c("navy", "white", "firebrick"))(100),
                   fontsize_row = 6,fontsize_col = 10)

df.anno <- data.frame(row.names = rownames(tmp))
df.anno$group <- tmp$dif_peak
table(df.anno$group)

pheatmap::pheatmap(zscore(log2(tmp[,-12]+1)), show_rownames = F, cluster_rows = F, cluster_cols = F,
                   col = colorRampPalette(c("navy", "white", "firebrick"))(100),
                   fontsize_row = 6,fontsize_col = 10, annotation_row = df.anno)


pheatmap::pheatmap(zscore(log2(tmp[,-c(12,1:2)]+1)), show_rownames = F, cluster_rows = F, cluster_cols = F,
                   col = colorRampPalette(c("navy", "white", "firebrick"))(100),
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
                   col = colorRampPalette(c("navy", "white", "firebrick"))(1000),
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

save_pheatmap_pdf(x, "heatmap.11.10.pdf", height = 30)


