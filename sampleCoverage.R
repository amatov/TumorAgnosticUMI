library("ggplot2")
library("gplots")
source("~/genomedk/matovanalysis/umiseq_analysis/R/read_bed.R")
pon_hg19 <- readRDS("~/genomedk/PolyA/faststorage/BACKUP/IMPROVE/call/references/200419_novaseq-xgen-sporacrc-pon.RDS") # 
str(pon_hg19)
# 45 Subjects of the Control Panel of Normal PON
pon_obj2 <- readRDS("sw_input_files/201020_hg38-novaseq-xgen-sporacrc-pon.RDS") # 
pon_counts <- pon_obj2[["pon"]]
no1 = array(0, dim=c(dim(pon_counts)[1],sum(list),dim(pon_counts)[3]))
for (i in 1:dim(pon_counts)[1]) {
  p2 <- data.frame(pon_counts[i,,])#PON
  p1 <- p2 [list == 1, ] 
  no1[i,,] <- data.matrix(p1)
}
no <- no1[,,1:4]+no1[,,6:9]
#########Gray band (variance) ############################
covP <- rowSums(no, dims = 2)
meaP<-apply(covP,2,mean) 
noP1 <- t(rbind(covP,meaP))
noP2 <- noP1[order(noP1[,(dim(pon_counts)[1]+1)],decreasing=TRUE),]
noP3 <- noP2[,1:dim(pon_counts)[1]]# the sorted coverage for each position in the 45 PONs
means = apply(noP3,1,mean)
sds = apply(noP3,1,sd)
df = data.frame(idx = 1:nrow(noP3), mean = means, sd = sds)
df2 = data.frame(sort(means, decreasing = TRUE)) 
ggplot(df, aes(x = idx, y = mean)) +geom_line() + geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd) , fill = "grey40",alpha=.5) 
 