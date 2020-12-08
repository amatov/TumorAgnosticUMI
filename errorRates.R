source("~/genomedk/matovanalysis/umiseq_analysis/R/read_bed.R") #1/0 list
f2 <- function(a, M, S){
  #Get indeces instead of number of sites
  which(apply(a, c(1, 2), function(x)sum(sum(x >= M, na.rm=T) >= S))>0)
}
#For one comination 
S =5; M=0.01; a = mafs
indP <- f2(mafs, M, S)
mafs[i] <- NA

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
# QIAGEN healthy samples ############################################################################################
pileupsQ <- list.files("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/qiagen_kit_test/201019", recursive = T, full.names = T, pattern = "bait.pileup")
countsQ0 <-  piles_to_counts(files = pileupsQ, 
                             regions = pon_hg19$regions)
countsQ= array(0, dim=c(dim(countsQ0)[1],sum(list),dim(countsQ0)[3]))
for (i in 1:dim(countsQ0)[1]) {
  p2 <- data.frame(countsQ0[i,,]) 
  p1 <- p2 [list == 1, ] 
  countsQ[i,,] <- data.matrix(p1)
}
countsQ1 <- countsQ[,,1:4] + countsQ[,,6:9]
Qc<-vector()
for (i in 1:24){
  Qc[i]<-sum(countsQ1[i,,])
}
erQ1<- vector()
mafsQ1 = array(0, dim=c(dim(countsQ1)[1],dim(countsQ1)[2],dim(countsQ1)[3]))
auxMQ <- rowSums(countsQ1, dims = 2) 
for (i in 1:dim(countsQ1)[1]) {
  #i=2
  mafsQ1[i,,] <- countsQ1[i,,]  /auxMQ[i,]
  erQ1[i] <- mean(mafsQ1[i,,][mafsQ1[i,,]<= 0.1])
}
plot(erQ1)
### IMPROVE #####
pileupsI <- list.files("~/genomedk/PolyA/faststorage/BACKUP/IMPROVE/sporacrc/N227", recursive = T, full.names = T, pattern = "bait.pileup")
countsI0 <-  piles_to_counts(files = pileupsI[1:184], 
                             regions = pon_hg19$regions)
countsI= array(0, dim=c(dim(countsI0)[1],sum(list),dim(countsI0)[3]))
for (i in 1:dim(countsI0)[1]) {
  p2 <- data.frame(countsI0[i,,])#CRUK
  p1 <- p2 [list == 1, ] 
  countsI[i,,] <- data.matrix(p1)
}
countsI1 <- countsI[,,1:4] + countsI[,,6:9]
erI1<- vector()
mafsI1 = array(0, dim=c(dim(countsI1)[1],dim(countsI1)[2],dim(countsI1)[3]))
auxMI <- rowSums(countsI1, dims = 2) 
for (i in 1:dim(countsI1)[1]) {
  mafsI1[i,,] <- countsI1[i,,]  /auxMI[i,]
  erI1[i] <- mean(mafsI1[i,,][mafsI1[i,,]<= 0.1])
}
plot(erI1)
# CRUK patient samples ##########################################################################################
pileupsC <- list.files("~/genomedk/PolyA/faststorage/BACKUP/CRUK/plasma/N289", recursive = T, full.names = T, pattern = "bait.pileup")
countsC0 <-  piles_to_counts(files = pileupsC[1:90], 
                             regions = pon_hg19$regions)
countsC= array(0, dim=c(dim(countsC0)[1],sum(list),dim(countsC0)[3]))
for (i in 1:dim(countsC0)[1]) {
  p2 <- data.frame(countsC0[i,,])#CRUK
  p1 <- p2 [list == 1, ] 
  countsC[i,,] <- data.matrix(p1)
}
countsC1 <- countsC[,,1:4] + countsC[,,6:9]
erC1<- vector()
mafsC1 = array(0, dim=c(dim(countsC1)[1],dim(countsC1)[2],dim(countsC1)[3]))
auxMC <- rowSums(countsC1, dims = 2) 
for (i in 1:dim(countsC1)[1]) {
  mafsC1[i,,] <- countsC1[i,,]  /auxMC[i,]
  erC1[i] <- mean(mafsC1[i,,][mafsC1[i,,]<= 0.1])
}
plot(erC1)
# DS samples ##################################################################################################
dat0 <- readRDS("sw_output_files/2020-10-23-145546_sw-output.RDS")  
pileupsD <- unlist(attributes(dat0))
pileupsD <- sub("/faststorage/project/PolyA/BACKUP", "~/genomedk/PolyA/faststorage/BACKUP", pileupsD)
all(file.exists(pileupsD))
countsD <-  piles_to_counts(files = pileupsD, 
                            regions = pon_obj2$regions)
countsDD= array(0, dim=c(dim(countsD)[1],sum(list),dim(countsD)[3]))
for (i in 1:dim(countsD)[1]) {
  p2 <- data.frame(countsD[i,,])
  p1 <- p2 [list == 1, ] 
  countsDD[i,,] <- data.matrix(p1)
}
countsD1 <- countsDD[,,1:4] + countsDD[,,6:9]
mafsD1 = array(0, dim=c(dim(countsD1)[1],dim(countsD1)[2],dim(countsD1)[3]))
erD1<- vector()
auxMD <- rowSums(countsD1, dims = 2) 
for (i in 1:dim(countsD1)[1]) {
  mafsD1[i,,] <- countsD1[i,,]  /auxMD[i,]
  erD1[i] <- mean(mafsD1[i,,][mafsD1[i,,]<= 0.1])
}
plot(erD1)
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