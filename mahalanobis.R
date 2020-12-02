
library(matrixStats)
library(reshape2)
#setwd ('G:\\PolyA/faststorage/BACKUP/N140_Targeting/specs/specs_analysis')
setwd ('~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/specs_analysis')
source("sw_input_files/duplex_tools.R")

library(ROCR)
source("U:\\Documents/R/utility_functions-master/recoder.R")
source("U:\\Documents/R/utility_functions-master/auc.R")
source("U:\\Documents/R/utility_functions-master/scaler.R")
source("U:\\Documents/R/utility_functions-master/confusion_plot.R")
source("G:\\PolyA/faststorage/BACKUP/N140_Targeting/specs/specs_analysis/sw_piles.R")
source("sw_input_files/duplex_tools.R")
library("dplyr")
library("tidyr")
source("~/genomedk/matovanalysis/umiseq_analysis/R/read_bed.R") #1/0 list
source("~/genomedk/matovanalysis/umiseq_analysis/R/mutationScore.R") #1/0 list
p2 <- data.frame(prior)#PON
p1 <- p2 [list == 1, ] 
prior2 <- data.matrix(p1)

# look up the 90 from this panel based on HG
# find the 15k positions for the 90
prior1 <- readRDS("sw_input_files/180903_prior.RDS")
prior <- prior1[,1:4]
# data frame to take list positions only
# source the updated R file mutationScore.R

pon_hg19 <- readRDS("~/genomedk/PolyA/faststorage/BACKUP/IMPROVE/call/references/200419_novaseq-xgen-sporacrc-pon.RDS") # 
str(pon_hg19)
# extract from COSMIC likelihood of the gene to be mutated and the mutation to be the same position as the mutation list
#cosmic_file <- read.table(file='~/genomedk/matovanalysis/umiseq_analysis/CosmicMutations/CosmicMutantExport.tsv',, sep = '\t', header = TRUE) #17G
cosmic_file <- read.table(file='~/genomedk/matovanalysis/umiseq_analysis/CosmicMutations/CosmicMutantExport_colon.tsv',, sep = '\t', header = TRUE) #0.6G dim(cosmic_file)1,511,817x40
#CosmicMutantExport_large_intestine # 1.7G

founderinfo <- readRDS("200330_founder-mutations-method1.RDS")
fI<-founderinfo[1:90,]# only the 90 mutations info
inMu1 <- fI$index# indexes of the 90 mutations on the panel
inMu <- vector()
inAu <- rep(0,18094*4)# empty vector 18094
inAu[inMu1]=1# put ones at the inMu positions
p2 <- data.frame(inAu) 
p1 <- p2 [list == 1, ] # we re loosing 4 of the 90 mutations based on the HG38 BED file and the HG19 mutation list
inMu <- which(p1>0)# multiply with "list"
list1 <- rep(0,15465*4)# empty vector 18094
list1[inMu]=1# put ones at the inMu positions
list2 <- matrix(list1, ncol = 4, byrow = 15465)
#muList <- p1
# 45 Subjects of the Control Panel of Normal PON ####################################################################
pon_obj2 <- readRDS("sw_input_files/201020_hg38-novaseq-xgen-sporacrc-pon.RDS") # 
pon_counts <- pon_obj2[["pon"]]
#str(pon_obj2)
#r1 <-rownames(pon_counts[1,,])
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
v<-apply(no,2:3,var) # variance of the counts of each position based on PON
#v<-apply(mafsP1,2:3,var) # variance of the VAFs of each position based on PON
v0 <- min(v[v>0])/10000000 # for counts w zero variance, we replace w a very small value
v1<- v
v1[v==0]=v0

#Test for normality the distributions of count values for each position over the 45 PON cohort
# if p<0.05 the distribution is not normal (ie its different from normal)
ponP = array(0, dim=c(dim(no)[2],dim(no)[3]))
sT = NULL
for (i in 1:dim(no)[2]){
  for (j in 1:dim(no)[3]){
    #i= 1 
    #j= 3
    #print(i)
    #print(j)
    if (sum(no[,i,j])==0){
      ponP[i,j]= 0
      }else{
    sT<- shapiro.test(no[,i,j])
    ponP[i,j] <- sT$p.value 
    }
  }
} #sum(ponP>=0.05)/61860: 3.5% of the positions are normally distributed over the PON cohort
# sum(ponP==0)/61860: 18.8% of the positions have always zeros for all subjects in the PON

# CRUK patient samples ##########################################################################################
pileupsC <- list.files("~/genomedk/PolyA/faststorage/BACKUP/CRUK/plasma/N289", recursive = T, full.names = T, pattern = "bait.pileup")
countsC0 <-  piles_to_counts(files = pileupsC[1:90], 
#countsC0 <-  piles_to_counts(files = pileupsC[1:8], 
                            regions = pon_hg19$regions)
countsC= array(0, dim=c(dim(countsC0)[1],sum(list),dim(countsC0)[3]))
for (i in 1:dim(countsC0)[1]) {
  #i=1
  p2 <- data.frame(countsC0[i,,])#CRUK
  p1 <- p2 [list == 1, ] 
  countsC[i,,] <- data.matrix(p1)
}
countsC1 <- countsC[,,1:4] + countsC[,,6:9]
covC <- rowSums(countsC1, dims = 2)
meaC<-apply(covC,2,mean) # variance of each position based on PON
plot(sort(meaC, decreasing = TRUE))
covC1 <- t(rbind(covC,meaC))
covC2 <- covC1[order(covC1[,91],decreasing=TRUE),]
plot(covC2[,91])# the sorted mean coverage for each position
erC1<- vector()
mafsC1 = array(0, dim=c(dim(countsC1)[1],dim(countsC1)[2],dim(countsC1)[3]))
auxMC <- rowSums(countsC1, dims = 2) 
for (i in 1:dim(countsC1)[1]) {
mafsC1[i,,] <- countsC1[i,,]  /auxMC[i,]
erC1[i] <- mean(mafsC1[i,,][mafsC1[i,,]<= 0.1])
}
plot(erC1)
mean(erC1)
w1C = array(0, dim=c(dim(countsC1)[1],dim(countsC1)[2],dim(countsC1)[3]))
w2C = array(0, dim=c(dim(countsC1)[1],dim(countsC1)[2],dim(countsC1)[3]))
sc1C <- vector()
sc2C <- vector()
for (i in 1:dim(mafsC1)[1]) {
  #i=1
  w1C[i,,] <- mafsC1[i,,]/v1/sum(list)  # the patient MAFs devided by the variance of the position based on PON'
  w2C[i,,] <- mafsC1[i,,]/v1*list2/sum(list2)  # the patient MAFs devided by the variance of the position based on PON'
  sc1C[i]<-sum(w1C[i,,])
  sc2C[i]<-sum(w2C[i,,])
}
plot(sc1C)
plot(sc2C)
#IMPROVE get the dates of OP
It <- read.table("~/genomedk/matovanalysis/umiseq_analysis/IMPROVEptList",header = TRUE)
table(It$op_time_cat)
tbI <- table(It$op_time_cat,It$pt_id)
tbI[1,] # preOP
tbI[2,] # postOP(14)
tbI[3,] # postOP(30)
It[,3] # IDs 
# CRUK samples #################################################################################################
pileupsC <- list.files("~/genomedk/PolyA/faststorage/BACKUP/CRUK/plasma/N289", recursive = T, full.names = T, pattern = "bait.pileup")
#countsIpre <-  piles_to_counts(files = pileupsC[9:90], 
                              countsIpre <-  piles_to_counts(files = pileupsC[1:8], 
                             regions = pon_hg19$regions)
# IMPROVE samples #################################################################################################
pileupsI <- list.files("~/genomedk/PolyA/faststorage/BACKUP/IMPROVE/sporacrc/N227", recursive = T, full.names = T, pattern = "bait.pileup")
It <- read.table("~/genomedk/matovanalysis/umiseq_analysis/IMPROVEptList",header = TRUE)
It$index <- sapply(as.character(It$library_id), function(x) grep(x, pileupsI)) # 141 of 179
preop_index <- It[ It$op_time_cat == -1, "index"] # 57
posop14_index <- It[ It$op_time_cat == 2, "index"] # 42
posop30_index <- It[ It$op_time_cat == 30, "index"] # 42

countsIpre <-  piles_to_counts(files = pileupsI[posop30_index], # only the first 178 samples are cancer
#countsI0 <-  piles_to_counts(files = pileupsI[1:184  ], # only the first 178 samples are cancer
                         regions = pon_hg19$regions)

#########function
#sc<-mahala(countsIpre, v1, list)  
########################
counts1= array(0, dim=c(dim(countsIpre)[1],sum(list),dim(countsIpre)[3]))
for (i in 1:dim(countsIpre)[1]) {
  p2 <- data.frame(countsIpre[i,,])#CRUK
  p1 <- p2 [list == 1, ] 
  counts1[i,,] <- data.matrix(p1)
}
counts <- counts1[,,1:4] + counts1[,,6:9]
mafsIpre= array(0, dim=c(dim(counts)[1],dim(counts)[2],dim(counts)[3]))
auxMIpre <- rowSums(counts, dims = 2) 
for (i in 1:dim(counts)[1]) {
  mafsIpre[i,,] <- counts[i,,]  /auxMIpre[i,]
}
w1Ipre = array(0, dim=c(dim(counts)[1],dim(counts)[2],dim(counts)[3]))
w11Ipre = array(0, dim=c(dim(counts)[1],dim(counts)[2],dim(counts)[3]))
w111Ipre = array(0, dim=c(dim(counts)[1],dim(counts)[2],dim(counts)[3]))
scIpre <- vector()
sc1Ipre <- vector()
sc11Ipre <- vector()
#mu<-apply(mafsP1,2:3,mean) # mean of the mafs of each position based on PON
for (i in 1:dim(mafsIpre)[1]) {
  #i=1
  w1Ipre[i,,] <- mafsIpre[i,,]/v1/(sum(list)*4)#prior[list]  # the patient MAFs devided by the variance of the position based on PON'
  scIpre[i] <- sum(w1Ipre[i,,])
  w11Ipre[i,,] <- mafsIpre[i,,][mafsIpre[i,,]<= 0.1]/v1/sum(mafsIpre[i,,]<= 0.1)  
  sc1Ipre[i] <- sum(w11Ipre[i,,])
  w111Ipre[i,,] <- mafsIpre[i,,][mafsIpre[i,,]<= 0.1]/v1/sum(mafsIpre[i,,]<= 0.1)*prior2
  #sc11Ipre[i]<-sum(w111Ipre[i,,])
  sc11Ipre[i]<-sum(sort(w111Ipre[i,,], decreasing = TRUE)[1:3])
}
scCctrl <- sc11Ipre # CR ctrl
scCpreOp <- sc11Ipre # CR preOp
scIpreOP <- sc11Ipre # IM preOP
scIpostOP14 <- sc11Ipre # IM postOP14
scIpostOP30 <- sc11Ipre # IM postOP30
#######################################################################
countsI0<-countsIpre
countsI= array(0, dim=c(dim(countsI0)[1],sum(list),dim(countsI0)[3]))
for (i in 1:dim(countsI0)[1]) {
  p2 <- data.frame(countsI0[i,,])#IMPROVE
  p1 <- p2 [list == 1, ] 
  countsI[i,,] <- data.matrix(p1)
}
countsI1 <- countsI[,,1:4] + countsI[,,6:9]
#covI <- rowSums(countsI1, dims = 2)
#meaI<-apply(covI,2,mean) # variance of each position based on PON
#plot(sort(meaI, decreasing = TRUE))
#covI1 <- t(rbind(covI,meaI))
#covI2 <- covI1[order(covI1[,180],decreasing=TRUE),]
#plot(covI2[,180])# the sorted mean coverage for each position
#mafsI1 <- countsI1/rowSums(countsI1)  
#mafsC1 <- countsC1/rowSums(countsC1)  
erI1<- vector()
mafsI1 = array(0, dim=c(dim(countsI1)[1],dim(countsI1)[2],dim(countsI1)[3]))
auxMI <- rowSums(countsI1, dims = 2) 
for (i in 1:dim(countsI1)[1]) {
  mafsI1[i,,] <- countsI1[i,,]  /auxMI[i,]
  erI1[i] <- mean(mafsI1[i,,][mafsI1[i,,]<= 0.1])
}
plot(erI1)
mean(erI1)

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
#mafsQ1 <- countsQ1/rowSums(countsQ1) 
erQ1<- vector()
mafsQ1 = array(0, dim=c(dim(countsQ1)[1],dim(countsQ1)[2],dim(countsQ1)[3]))
auxMQ <- rowSums(countsQ1, dims = 2) 
for (i in 1:dim(countsQ1)[1]) {
  #i=2
  mafsQ1[i,,] <- countsQ1[i,,]  /auxMQ[i,]
  erQ1[i] <- mean(mafsQ1[i,,][mafsQ1[i,,]<= 0.1])
}
plot(erQ1)
mean(erQ1)
w1Q = array(0, dim=c(dim(countsQ1)[1],dim(countsQ1)[2],dim(countsQ1)[3]))
w11Q = array(0, dim=c(dim(countsQ1)[1],dim(countsQ1)[2],dim(countsQ1)[3]))
w111Q = array(0, dim=c(dim(countsQ1)[1],dim(countsQ1)[2],dim(countsQ1)[3]))
w2Q = array(0, dim=c(dim(countsQ1)[1],dim(countsQ1)[2],dim(countsQ1)[3]))
sc1Q <- vector()
sc11Q <- vector()
sc111Q <- vector()
sc2Q <- vector()
for (i in 1:dim(mafsQ1)[1]) {
  #i=21
  w1Q[i,,] <- mafsQ1[i,,]/v1/(sum(list)*4)  # the patient MAFs devided by the variance of the position based on PON'
  sc1Q[i]<-sum(w1Q[i,,])
  #w2Q[i,,] <- mafsQ1[i,,]/v1*list2/sum(list2) # the patient MAFs devided by the variance of the position based on PON'
  sc2Q[i]<-sum(w2Q[i,,])
  w11Q[i,,] <- mafsQ1[i,,][mafsQ1[i,,]<= 0.1]/v1/sum(mafsQ1[i,,]<= 0.1) 
  sc11Q[i]<-sum(w11Q[i,,])
  w111Q[i,,] <- mafsQ1[i,,][mafsQ1[i,,]<= 0.1]/v1/sum(mafsQ1[i,,]<= 0.1)*prior2
  #sc111Q[i]<-sum(w111Q[i,,])
  sc111Q[i]<-sum(sort(w111Q[i,,], decreasing = TRUE)[1:3])
}
r1<- -3
r2 <- 14
plot(log2(sc1Q), ylim=range(c(r1,r2)), col="green", pch = 17)
par(new = TRUE)
plot(log2(sc11Q), ylim=range(c(r1,r2)), col="red", pch = 19)
par(new = TRUE)
plot(log2(sc111Q), ylim=range(c(r1,r2)), col="blue", pch = 18)
#plot(sc2Q)
# DS samples ##################################################################################################
dat0 <- readRDS("sw_output_files/2020-10-23-145546_sw-output.RDS")  
pileupsD <- unlist(attributes(dat0))
pileupsD <- sub("/faststorage/project/PolyA/BACKUP", "~/genomedk/PolyA/faststorage/BACKUP", pileupsD)
all(file.exists(pileupsD))
countsD <-  piles_to_counts(files = pileupsD, 
                            regions = pon_obj2$regions)
countsDD= array(0, dim=c(dim(countsD)[1],sum(list),dim(countsD)[3]))
for (i in 1:dim(countsD)[1]) {
  p2 <- data.frame(countsD[i,,])#IMPROVE
  p1 <- p2 [list == 1, ] 
  countsDD[i,,] <- data.matrix(p1)
}
countsD1 <- countsDD[,,1:4] + countsDD[,,6:9]
erD1<- vector()
auxMD <- rowSums(countsD1, dims = 2) 
for (i in 1:dim(countsD1)[1]) {
  mafsD1[i,,] <- countsD1[i,,]  /auxMD[i,]
  erD1[i] <- mean(mafsD1[i,,][mafsD1[i,,]<= 0.1])
}
plot(erD1)
mean(erD1)
w1D = array(0, dim=c(dim(countsD1)[1],dim(countsD1)[2],dim(countsD1)[3]))
w2D = array(0, dim=c(dim(countsD1)[1],dim(countsD1)[2],dim(countsD1)[3]))
sc1D <- vector()
sc2D <- vector()
for (i in 1:dim(mafsD1)[1]) {
  #i=1
  w1D[i,,] <- mafsD1[i,,]/v1/sum(list)#*10000 # the patient MAFs devided by the variance of the position based on PON'
  sc1D[i]<-sum(w1D[i,,])
  w2D[i,,] <- mafsD1[i,,]/v1*list2/sum(list2)  # the patient MAFs devided by the variance of the position based on PON'
  sc2D[i]<-sum(w2D[i,,])
}
plot(sc1D)
plot(sc2D)
#################tumor-informed score########################################################################
r1<- -27
r2 <- -3
plot(log2(sc2I), ylim=range(c(r1,r2)), col="orange", main = "Tumor-Informed Score", pch = 15)
par(new = TRUE)
plot(log2(sc2C), ylim=range(c(r1,r2)), col="red", pch = 19)
par(new = TRUE)
plot(log2(sc2Q), ylim=range(c(r1,r2)), col="green", pch = 17)
par(new = TRUE)
plot(log2(sc2D), ylim=range(c(r1,r2)), col="blue", pch = 11)
legend(3.05,3.5,legend=c("IMPROVE", "CRUK","QIAGEN", "DS"),col=c("orange","red","green","blue"),lty=1:1, cex=1.0)
##################tumor non-informed score#########################################################################################
r1<- -8
r2 <- -1
plot(log2(scIpreOP), ylim=range(c(r1,r2)), col="red", main = "Score", pch = 15)
par(new = TRUE)
plot(log2(scCpreOp), ylim=range(c(r1,r2)), col="red", pch = 19)
par(new = TRUE)
plot(log2(scCctrl), ylim=range(c(r1,r2)), col="green", pch = 19)
par(new = TRUE)
plot(log2(scIpostOP14), ylim=range(c(r1,r2)), col="green", pch = 15)
par(new = TRUE)
plot(log2(scIpostOP30), ylim=range(c(r1,r2)), col="yellow", pch = 15)
par(new = TRUE)
plot(log2(sc111Q), ylim=range(c(r1,r2)), col="blue", pch = 15)
legend(18,-0.1,legend=c("IMpreOp sq", "CRpreOp ci","CRctrl ci","IM14 sq", "IM30 sq", "Qia ci"),col=c("red","red","green","green","yellow","blue"),lty=1:1, cex=1.0)
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
###
plot(log2(erI1), ylim=range(c(-15.9,-13.1)), col="orange", main = "Error rate", pch = 15)
par(new = TRUE)
plot(log2(erC1), ylim=range(c(-15.9,-13.1)), col="red", pch = 19)
par(new = TRUE)
plot(log2(erQ1), ylim=range(c(-15.9,-13.1)), col="green", pch = 17)
par(new = TRUE)
plot(log2(erP1), ylim=range(c(-15.9,-13.1)), col="purple", pch = 12)
par(new = TRUE)
plot(log2(erD1), ylim=range(c(-15.9,-13.1)), col="blue", pch = 11)
legend(3,1.1e-04,legend=c("IMPROVE", "CRUK","QIAGEN", "PON","DS"),col=c("orange","red","green","purple","blue"),lty=1:1, cex=1.0)
# merge sc1C and sc1I with different tags #################################################################
df1 <- data.frame(sc1I)# IMPROVE 178
df1$cohort = 1
df2 <- data.frame(sc1C)# CRUK 82
df2$cohort = 2
df <- bind_rows(df1, df2) # change names to make sure they merge
# sort and plot

mahala <- function(counts0, v, list) {
  
  counts1= array(0, dim=c(dim(counts0)[1],sum(list),dim(counts0)[3]))
  for (i in 1:dim(counts)[1]) {
    p2 <- data.frame(counts0[i,,])#CRUK
    p1 <- p2 [list == 1, ] 
    counts1[i,,] <- data.matrix(p1)
  }
  counts <- counts1[,,1:4] + counts1[,,6:9]
  
  mafs = array(0, dim=c(dim(counts)[1],dim(counts)[2],dim(counts)[3]))
  auxM <- rowSums(counts, dims = 2) 
  for (i in 1:dim(counts)[1]) {
    mafs[i,,] <- counts[i,,]  /auxM[i,]
  }
  w1 = array(0, dim=c(dim(counts)[1],dim(counts)[2],dim(counts)[3]))
  sc <- vector()
  for (i in 1:dim(mafs)[1]) {
    print(i)
    w1[i,,] <- mafs[i,,]#/v/sum(list)  # the patient MAFs devided by the variance of the position based on PON'
    sc[i]<-sum(w1[i,,])
  }
  return(sc) # mafs, counts
}
