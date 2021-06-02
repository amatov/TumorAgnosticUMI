library("dplyr")
library("glmnet")
source("~/genomedk/matovanalysis/umiseq_analysis/R/read_bed.R") #1/0 list
library("ggplot2")
library("precrec")
library("ROCit")
# add new flags; number of fragments, age, gender, concentration, 10 flags
library(GenomicRanges)
source("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/specs_analysis/sw_input_files/tools.R")
#################################################
w3_files <- list.files("~/genomedk/matovanalysis/umiseq_analysis/CRUK5Mb", recursive = T, full.names = T, pattern = "_consensus.txt")
W3_N289_70 <- read.table("~/genomedk/matovanalysis/umiseq_analysis/CRUK5Mb/N289-70.txt", header = T)  
W3_N289_85 <- read.table("~/genomedk/matovanalysis/umiseq_analysis/CRUK5Mb/C44A06969D_cfdna_N289_85_consensus.txt")#, header = T)  
W3_N289_58 <- read.table("~/genomedk/matovanalysis/umiseq_analysis/CRUK5Mb/C47A07007D_cfdna_N289_58_consensus.txt")#, header = T)  
w3 <- W3_N289_58[,3:702] #mean(unlist(w3)) 0.007291663

w3 <- W3_N289_70[,3:702]
W3_N289_85
w3 <- W3_N289_85[,3:702]

mean(unlist(w3))# 0.007291663
w31<-unlist(w3)
w32<-w31[which(w31> 0)]
muPos <- sum(w31) / length(w32) # 332.431

#w3t <- unlist(w3[219,])
#w3t1 <- unlist(w3t[w3t>0])
length(w3t1) #232
muPos <- sum(w3t) / length(w3t1) # 332.431
meanS <- mean(muPos)
meanS # 7.065037

muPos <- vector()
w3pos <- vector()
w3pos2 <- vector()

for (i in 1:length(w32)){
  #i=219 # 219, 220, 221, 222
  w3t <- unlist(w3[i,])
  w3t1 <- w3t[which(w3t>0)]

  length(w3t1) 
  auxPos <- (meanS - w3t1)*(meanS - w3t1)
  w3pos[i] <- mean(auxPos)
  muPos[i] <- sum(w3t) / length(w3t1) # 332.431
  auxPos2 <- (muPos[i] - w3t1)*(muPos[i] - w3t1)
  
  w3pos2[i] <- mean(auxPos2)
  
}
############################################################################
countsW <- readRDS("~/genomedk/matovanalysis/umiseq_analysis/R/cruki.RDS") # 

# 45 Subjects of the Control Panel of Normal PON ####################################################################
pon_obj2 <- readRDS("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/umiseq_paper/reference/201217_hg38-novaseq-xgen-sporacrc-pon.RDS") # 46
pon_counts <- pon_obj2[["pon"]]
no0 = array(0, dim=c(dim(pon_counts)[1]-1,dim(pon_counts)[2],dim(pon_counts)[3]))
no0[1:27,,] <- pon_counts[1:27,,]
no0[28:45,,]<-pon_counts[29:46,,]

no <- no0[,,1:4]+no0[,,6:9]
mafsP1 = array(0, dim=c(dim(no)[1],dim(no)[2],dim(no)[3]))
auxMP <- rowSums(no, dims = 2) 
for (i in 1:dim(no)[1]) {
  mafsP1[i,,] <- no[i,,]  /auxMP[i,]
}
############################# Standard Variation of the VAFs of PON at each of the 62k positions #######################################3
#v<-apply(mafsP1,2:3,sd) # variance of the counts of each position based on PON
#v0 <- min(v[v>0])/100000#00 # for counts w zero variance, we replace w a very small value
#v1<- v
#v1[v==0]=v0
v1 <- apply(mafsP1,2:3,sd, na.rm = T) #Note na.rm = T
v1[ v1 == 0 ] <- 1E-5
# PON mutations and variability ###################################################################################
sitemut <- t(apply(pon_obj2$coordinates, 1, function(x){
  paste( paste0(trimws(x[1]), ":", trimws(x[2]), "_"),
         paste(x[3] ,c("A", "T", "C", "G"), sep = "/"),
         sep = "")}))
########################
prior0 <- readRDS("~/genomedk/matovanalysis/umiseq_analysis/R/cosmic-prior.RDS") # 

# cruk plasma data ##################################################################################################
countsC00 <- readRDS("~/genomedk/matovanalysis/umiseq_analysis/R/cruk-counts.RDS") # 
countsC001 <- countsC00[,,1:4] + countsC00[,,6:9]
mafsC = array(0, dim=c(dim(countsC001)[1],dim(countsC001)[2],dim(countsC001)[3]))
auxC <- rowSums(countsC001, dims = 2) 
for (i in 1:dim(countsC001)[1]) {
  mafsC[i,,] <- countsC001[i,,]  /auxC[i,]#/v1#/sum(list) # Weight 1
}

countsC0= array(0, dim=c(95,dim(countsC00)[2],dim(countsC00)[3]))
countsC0[1:69,,]<-countsC00[2:70,,] # take only cancer samples 2:70
countsC0[70:95,,]<-countsC00[105:130,,] # take only cancer samples 105:130
countsC= array(0, dim=c(dim(countsC0)[1],sum(list),dim(countsC0)[3]))

countsC1 <- countsC[,,1:4] + countsC[,,6:9]
mafsC1 = array(0, dim=c(dim(countsC1)[1],sum(list),dim(countsC1)[3]))
mafsC2= array(0, dim=c(dim(countsC1)[1],sum(list)*dim(countsC1)[3]))
auxMC <- rowSums(countsC1, dims = 2) 
for (i in 1:dim(countsC1)[1]) {
  mafsC1[i,,] <- countsC1[i,,]  /auxMC[i,]#/v1/sum(list) # Weight 1
}

cancer_list <- unique(cancer_SNPs$i)
noCancer_list <- c(adenoma_list, control_list)

mahaT <- vector()
j=1
for (i in cancer_list){
  for (k in control_list[6:15]){
    mahaT[j] <- sum(countsW[i,,]* mafsC[k,,]/v1,  na.rm=T)/sum(countsW[i,,])#/sum(mafsC[k,,],  na.rm=T)
    print(sum(mafsC[k,,],  na.rm=T))
    j = j + 1
  }
}
length(mahaT)
plot(mahaT)

mahaT2 <- vector()
j=1
for (i in cancer_list){
  for (k in adenoma_list){
    mahaT2[j] <- sum(countsW[i,,]* mafsC[k,,]/v1,  na.rm=T)/sum(countsW[i,,])#/sum(mafsC[k,,],  na.rm=T) # Weight 2
    j = j + 1
  }
}
length(mahaT2)
plot(mahaT2)

maha <- vector()
j=1
for (i in cancer_list){
  maha[j] <- sum(countsW[i,,]* mafsC[i,,]/v1,  na.rm=T) / sum(countsW[i,,])#/sum(mafsC[i,,],  na.rm=T)#45 of 95 cancers have WES, of them 37 have VAF>0

if (maha[j]>0) {
j = j + 1
  }

mahaCancer <- maha[1:40] # 37 of 45 , w 8 cancers w VAFs 0
plot(mahaCancer) # Cancer samples with VAF>0
mahaControl <- c(mahaT)#, mahaT2) # 1575 (45x15+45x20)
plot(mahaT2,ylim=range(c(0,0.0035)))
plot(mahaT,ylim=range(c(0,0.0035)))
#ROC; find where overall success rate numbers are
rocTF = rep(0, (40+45*10))#15))#+45*20))
rocTF[1:40]=1
##################################
pred <- prediction(c(mahaCancer, mahaControl), rocTF)
perf<-performance(pred,"tpr", "fpr")
plot(perf)
auc<- performance(pred,"auc")
auc
##########################AUC#######################################################################
pROC_obj <- roc(rocTF,c(mahaCancer, mahaControl),
                smoothed = TRUE,
                # arguments for ci
                ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                # arguments for plot
                plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                print.auc=TRUE, show.thres=TRUE)

sens.ci <- ci.se(pROC_obj)
plot(sens.ci, type="shape", col="lightblue")
## Warning in plot.ci.se(sens.ci, type = "shape", col = "lightblue"): Low
## definition shape.
plot(sens.ci, type="bars")
################################################################
ROCit_obj <- rocit(score=c(mahaCancer, mahaControl),class=rocTF)
plot(ROCit_obj)
ksplot(ROCit_obj)
###############
precrec_obj <- evalmod(scores = c(mahaCancer, mahaControl), labels = rocTF)
autoplot(precrec_obj)
precrec_obj2 <- evalmod(scores = c(mahaCancer, mahaControl), labels = rocTF, mode="basic")
autoplot(precrec_obj2)  


############## RDS 95 cancer, 20 adenoma, 15 control #####################################################
sitemutPON # mutation map 15k x 4 for each position of the core panel
length(plasma$sitemut_hg38) # 294 WES mutations
#plasma$i <- sapply(as.character(plasma$sitemut_hg38), function(s) which(sitemutPON %in% s))

# loop over 95 cancers (mafsC1)95,15k,4 
cancer_SNPs$sample_label # list of samples with WES data
cancer_SNPs$sitemut_hg38 # list of WES mutations
cancerSamples <- unique(cancer_SNPs$sample_label)

dimnames(countsC00)[[1]] #The pileup names are the first dim
indexInRDS <- sapply(cancerSamples, function(x) grep(x, pileupsWES))  
indexInRDS[2] #C22A06329D  6 
cancerSamples[1] # "C22A06329D" There are 45 cancer CRUK samples with WES

scores <- array(0,dim=c(dim(mafsC1)[1],dim(mafsC1)[1]))
# test for cancer1 the score for its 3 mutations compared to all 95 cancers at the same mutations

plasma$ind <- sapply(as.character(sitemutPON), function(s) which(cancer_SNPs$sitemut_hg38 %in% s)) # INDEX IN COUNTS MATRIX OF MUTATIONS
plasma$i <- sapply(as.character(plasma$sitemut_hg38), function(s) which(sitemut %in% s))



