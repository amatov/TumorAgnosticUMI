source("~/genomedk/matovanalysis/umiseq_analysis/R/read_bed.R") #1/0 list
f2 <- function(a, M, S){
  #Get indeces instead of number of sites
  which(apply(a, c(1, 2), function(x)sum(sum(x >= M, na.rm=T) >= S))>0)
}

#For one comination 
S =5; M=0.01; a = mafs
indP <- f2(mafs, M, S)
# 45 Subjects of the Control Panel of Normal PON ####################################################################
pon_obj2 <- readRDS("sw_input_files/201020_hg38-novaseq-xgen-sporacrc-pon.RDS") # 
pon_counts <- pon_obj2[["pon"]]
no1 = array(0, dim=c(dim(pon_counts)[1],sum(list),dim(pon_counts)[3]))
for (i in 1:dim(pon_counts)[1]) {
  p2 <- data.frame(pon_counts[i,,])#PON
  p1 <- p2 [list == 1, ] 
  no1[i,,] <- data.matrix(p1)
}
no <- no1[,,1:4]+no1[,,6:9]
mafsP1 = array(0, dim=c(dim(no)[1],dim(no)[2],dim(no)[3]))
erP1<- vector()
auxMP <- rowSums(no, dims = 2) 
for (i in 1:dim(no)[1]) {
  #i = 1
  mafsP1[i,,] <- no[i,,]  /auxMP[i,]
  erP1[i] <- mean(mafsP1[i,,][mafsP1[i,,]<= 0.1])
}
plot(erP1)
mean(erP1)















###########################################################################################################
plot(erI1, ylim=range(c(1e-05,0.00011)), col="orange", main = "Error rate", pch = 15)
par(new = TRUE)
plot(erC1, ylim=range(c(1e-05,0.00011)), col="red", pch = 19)
par(new = TRUE)
plot(erQ1, ylim=range(c(1e-05,0.00011)), col="green", pch = 17)
par(new = TRUE)
plot(erP1, ylim=range(c(1e-05,0.00011)), col="purple", pch = 12)
par(new = TRUE)
plot(erD1, ylim=range(c(1e-05,0.00011)), col="blue", pch = 11)
legend(3,1.1e-04,legend=c("IMPROVE", "CRUK","QIAGEN", "PON","DS"),col=c("orange","red","green","purple","blue"),lty=1:1, cex=1.0)