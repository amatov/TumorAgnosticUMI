library("dplyr")
library("glmnet")
source("~/genomedk/matovanalysis/umiseq_analysis/R/read_bed.R") #1/0 list
library("ggplot2")
library("precrec")
library("ROCit")
# add new flags; number of fragments, age, gender, concentration, 10 flags
library(GenomicRanges)
source("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/specs_analysis/sw_input_files/tools.R")
# 45 Subjects of the Control Panel of Normal PON ####################################################################
pon_obj2 <- readRDS("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/umiseq_paper/reference/201217_hg38-novaseq-xgen-sporacrc-pon.RDS") # 46
pon_counts <- pon_obj2[["pon"]]
no0 = array(0, dim=c(dim(pon_counts)[1]-1,dim(pon_counts)[2],dim(pon_counts)[3]))
no0[1:27,,] <- pon_counts[1:27,,]
no0[28:45,,]<-pon_counts[29:46,,]
#no1 = array(0, dim=c(dim(no0)[1],sum(list),dim(no0)[3]))
#for (i in 1:dim(no0)[1]) {
  #i=1
#  p2 <- data.frame(no0[i,,])#PON
#  p1 <- p2 [list == 1, ] 
#  no1[i,,] <- data.matrix(p1)
#}
no <- no0[,,1:4]+no0[,,6:9]
mafsP1 = array(0, dim=c(dim(no)[1],dim(no)[2],dim(no)[3]))
auxMP <- rowSums(no, dims = 2) 
for (i in 1:dim(no)[1]) {
  mafsP1[i,,] <- no[i,,]  /auxMP[i,]
}
############################# Variance of the PON at each of the 62k positions #######################################3
v<-apply(mafsP1,2:3,sd) # variance of the counts of each position based on PON
v0 <- min(v[v>0])/100000#00 # for counts w zero variance, we replace w a very small value
v1<- v
v1[v==0]=v0
# PON mutations and variability ###################################################################################
sitemut <- t(apply(pon_obj2$coordinates, 1, function(x){
  paste( paste0(trimws(x[1]), ":", trimws(x[2]), "_"),
         paste(x[3] ,c("A", "T", "C", "G"), sep = "/"),
         sep = "")}))
#dfS <- data.frame(sitemut) 
#sitemutPON <- dfS[list == 1,] 
###########QIAGEN
# QIAGEN healthy samples ############################################################################################
pileupsQ <- list.files("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/qiagen_kit_test/201019", recursive = T, full.names = T, pattern = "bait.pileup")
countsQ00 <-  piles_to_counts(files = pileupsQ, 
                              regions = pon_obj2$coordinates)
#countsQ = array(0, dim=c(dim(countsQ00)[1]-2,dim(countsQ00)[2],dim(countsQ00)[3]))
countsQ = array(0, dim=c(dim(countsQ00)[1]-1,dim(countsQ00)[2],dim(countsQ00)[3]))
countsQ[1:16,,] <- countsQ00[1:16,,]
#countsQ[17:22,,]<-countsQ00[19:24,,]
countsQ[17:23,,]<-countsQ00[18:24,,]
countsQ0 <- countsQ
countsQ= array(0, dim=c(dim(countsQ0)[1],sum(list),dim(countsQ0)[3]))
for (i in 1:dim(countsQ0)[1]) {
  countsQ0[i,,][indP]<-NA # blacklist
  p2 <- data.frame(countsQ0[i,,]) 
  p1 <- p2 [list == 1, ] 
  countsQ[i,,] <- data.matrix(p1)
}
countsQ1 <- countsQ[,,1:4] + countsQ[,,6:9]  
erQ1<- vector()
mafsQ1 = array(0, dim=c(dim(countsQ1)[1],dim(countsQ1)[2],dim(countsQ1)[3]))
auxMQ <- rowSums(countsQ1, dims = 2) 
for (i in 1:dim(countsQ1)[1]) {
  #i=2
  mafsQ1[i,,] <- countsQ1[i,,]  /auxMQ[i,]
  erQ1[i] <- mean(mafsQ1[i,,][mafsQ1[i,,]<= VAFcut], na.rm=T)  
}
plot(erQ1)
print(erQ1)
# cruk plasma data ##################################################################################################
countsC00 <- readRDS("~/genomedk/matovanalysis/umiseq_analysis/R/cruk-counts.RDS") # 
countsC001 <- countsC00[,,1:4] + countsC00[,,6:9]
mafsC = array(0, dim=c(dim(countsC001)[1],dim(countsC001)[2],dim(countsC001)[3]))
auxC <- rowSums(countsC001, dims = 2) 
for (i in 1:dim(countsC001)[1]) {
  mafsC[i,,] <- countsC001[i,,]  /auxC[i,]#/v1#/sum(list) # Weight 1
  #mafsC2[i,] <- mafsC1[i,,]/sum(mafsC1[i,,]) # Normalize 
}
dimnames(countsC00)[[1]] #The pileup names are the first dim
#listSamples <- unlist(sapply(cancer_SNPs$sample_label , function(x) grep(x, x = dimnames(countsC00)[1] )))
pileupsWES <- dimnames(countsC00)[[1]][1:95] 
pileupsWES [1:69]<- dimnames(countsC00)[[1]][2:70] 
pileupsWES [70:95] <- dimnames(countsC00)[[1]][105:130]

countsC0= array(0, dim=c(95,dim(countsC00)[2],dim(countsC00)[3]))
countsC0[1:69,,]<-countsC00[2:70,,] # take only cancer samples 2:70
countsC0[70:95,,]<-countsC00[105:130,,] # take only cancer samples 105:130
countsC= array(0, dim=c(dim(countsC0)[1],sum(list),dim(countsC0)[3]))
#for (i in 1:dim(countsC0)[1]) {
  #i=1
#  p2 <- data.frame(countsC0[i,,])#CRUK
#  p1 <- p2 [list == 1, ] 
#  countsC[i,,] <- data.matrix(p1)
#}
countsC1 <- countsC[,,1:4] + countsC[,,6:9]
mafsC1 = array(0, dim=c(dim(countsC1)[1],sum(list),dim(countsC1)[3]))
mafsC2= array(0, dim=c(dim(countsC1)[1],sum(list)*dim(countsC1)[3]))
auxMC <- rowSums(countsC1, dims = 2) 
for (i in 1:dim(countsC1)[1]) {
  mafsC1[i,,] <- countsC1[i,,]  /auxMC[i,]#/v1/sum(list) # Weight 1
  #mafsC2[i,] <- mafsC1[i,,]/sum(mafsC1[i,,]) # Normalize 
}

#C22A06329D - #6, 3 mutations
#C32A06485D - #10, 4 mutations
cancer_SNPs$i # 45 unique
#  6  6  6  7  8  9  9 10 10 10 10 11 11 11 12 12 12 12 18 18 18 18 18 21 21 21 22 22 22 22 23 23 24 24 24 24 24 25 25 26 26 26 26 27 27
# 27 28 28 31 31 31 32 32 32 32 33 33 33 35 35 37 37 37 37 38 38 39 39 41 43 43 43 44 44 44 45 46 46 48 48 49 49 50 50 50 51 51 54 54 54
# 54 55 55 56 56 58 58 58 60 60 60 62 63 64 65 66 66 67 67 67 68 68 68 69
cancer_list <- unique(cancer_SNPs$i)
noCancer_list <- c(adenoma_list, control_list)

mahaT <- vector()
j=1
for (i in cancer_list){
  #i = 6
  for (k in control_list){
    mahaT[j] <- sum(countsW[i,,]* mafsC[k,,]/v1,  na.rm=T)/sum(countsW[i,,])/sum(mafsC[k,,],  na.rm=T)
    j = j + 1
  }
  #print(mahaT)
}
length(mahaT)
plot(mahaT)

mahaT2 <- vector()
j=1
for (i in cancer_list){
  #i = 6
  for (k in adenoma_list){
    mahaT2[j] <- sum(countsW[i,,]* mafsC[k,,]/v1,  na.rm=T)/sum(countsW[i,,])/sum(mafsC[k,,],  na.rm=T)
    j = j + 1
  }
  #print(mahaT2)
}
length(mahaT2)
plot(mahaT2)

maha <- vector()
j=1
for (i in cancer_list){
  maha[j] <- sum(countsW[i,,]* mafsC[i,,]/v1,  na.rm=T) / sum(countsW[i,,])/sum(mafsC[k,,],  na.rm=T)#45 of 95 cancers have WES, of them 37 have VAF>0
  print(sum(countsW[i,,]))
# maha[i] <- sum(countsW[6,,]* mafsC[i,,]/v1,  na.rm=T) # should be VAF mafsC[i,,]) #
if (maha[j]>0) {
j = j + 1
  }
#aux <- sum(countsW[6,,]*countsC001[i,,]/v1,  na.rm=T) # should be VAF mafsC[i,,]) #
# for countsC001 #6 [1] 142
# devide by PON var each pos
} 
print(maha)

mahaCancer <- maha[1:37] # 37 of 45 , w 8 cancers w VAFs 0
plot(mahaCancer) # Cancer samples with VAF>0

mahaControl <- c(mahaT)#, mahaT2) # 1575 (45x15+45x20)
plot(mahaT2,ylim=range(c(0,0.0035)))
plot(mahaT,ylim=range(c(0,0.0035)))
#ROC; find where overall success rate numbers are
rocTF = rep(0, (37+45*15))#+45*20))
rocTF[1:37]=1

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
#######################
# highest Youden’s score index 9: max (sensitivity + specificity-1) = 1.6568 for .9 specificity/1-FP(10%) and .7568 sensitivity/TP(75.68%)
sort(mahaT,decreasing = T)[1:68]
#[1] 1.143055e-03 9.922270e-04 5.618856e-04 4.961135e-04 4.236979e-04 3.307423e-04 3.236866e-04 2.979202e-04 2.749379e-04 2.650314e-04 2.606588e-04
#[12] 2.604915e-04 2.596182e-04 2.468367e-04 2.176489e-04 2.162913e-04 2.052583e-04 1.822990e-04 1.811689e-04 1.778197e-04 1.746562e-04 1.726651e-04
#[23] 1.700331e-04 1.632367e-04 1.606742e-04 1.591930e-04 1.580977e-04 1.575442e-04 1.566605e-04 1.529997e-04 1.489601e-04 1.489601e-04 1.433361e-04
#[34] 1.429947e-04 1.419545e-04 1.417240e-04 1.415791e-04 1.321530e-04 1.298091e-04 1.252298e-04 1.221734e-04 1.217086e-04 1.210090e-04 1.206197e-04
#[45] 1.198209e-04 1.180883e-04 1.160795e-04 1.133554e-04 1.116930e-04 1.116671e-04 1.065573e-04 1.056450e-04 1.054419e-04 1.033838e-04 1.021434e-04
#[56] 1.014164e-04 9.896396e-05 9.889234e-05 9.798165e-05 9.604065e-05 9.542160e-05 9.463632e-05 9.463632e-05 9.463632e-05 9.325567e-05 9.298205e-05
#[67] 8.833470e-05 8.833470e-05
sort(mahaCancer)[1:30]
#[1] 2.039784e-05 2.511039e-05 2.782754e-05 2.843320e-05 2.871902e-05 3.929044e-05 4.701047e-05 5.759859e-05 7.149662e-05 1.057386e-04 1.187934e-04
#[12] 1.196277e-04 1.274658e-04 1.614950e-04 1.854260e-04 1.854770e-04 2.166700e-04 2.728979e-04 2.984374e-04 3.272750e-04 3.321927e-04 3.453298e-04
#[23] 3.588652e-04 4.054661e-04 4.082032e-04 4.092563e-04 4.339736e-04 5.271311e-04

#For a threshold: .72 specificity/1-FP(28%, 189 of 675) for controls + .973 sensitivity/TP(97.3%, 36 of 37) for cancers. 
#95% CI (2000 stratified bootstrap replicates): 
ci.thresholds(pROC_obj)$specificity
ci.thresholds(pROC_obj)$sensitivity

dimnames(countsC00)[[1]][cancer_list[which(mahaFull==0)]]
#[1] "/faststorage/project/PolyA/BACKUP/CRUK/plasma/N289/6678/C28A06678D_cfdna_N289-21___200630bix195220/output/C28A06678D_cfdna_N289-21_consensus.bait.pileup"
#[2] "/faststorage/project/PolyA/BACKUP/CRUK/plasma/N289/6721/C23A06721D_cfdna_N289-25___200630vyr195214/output/C23A06721D_cfdna_N289-25_consensus.bait.pileup"
#[3] "/faststorage/project/PolyA/BACKUP/CRUK/plasma/N289/6780/C36A06780D_cfdna_N289-31___200630myl195228/output/C36A06780D_cfdna_N289-31_consensus.bait.pileup"
#[4] "/faststorage/project/PolyA/BACKUP/CRUK/plasma/N289/6792/C62A06792D_cfdna_N289-33___200630qev195246/output/C62A06792D_cfdna_N289-33_consensus.bait.pileup"
#[5] "/faststorage/project/PolyA/BACKUP/CRUK/plasma/N289/6901/C54A06901D_cfdna_N289-41___200630wuj195243/output/C54A06901D_cfdna_N289-41_consensus.bait.pileup"
#[6] "/faststorage/project/PolyA/BACKUP/CRUK/plasma/N289/7037/C43A07037D_cfdna_N289-63___200630wil195235/output/C43A07037D_cfdna_N289-63_consensus.bait.pileup"
#[7] "/faststorage/project/PolyA/BACKUP/CRUK/plasma/N289/7117/C17A07117D_cfdna_N289-71___200630pit195212/output/C17A07117D_cfdna_N289-71_consensus.bait.pileup"
#[8] "/faststorage/project/PolyA/BACKUP/CRUK/plasma/N289/7334/C42A07334D_cfdna_N289-87___200921zan150711/output/C42A07334D_cfdna_N289-87_consensus.bait.pileup"

# control#1 is also HG19 
#[1] "/faststorage/project/PolyA/BACKUP/CRUK/plasma/N289/421486/S07A05776D_cfdna_N289-101___200921wyw150716/output/S07A05776D_cfdna_N289-101_consensus.bait.pileup"

##############################SW##################################################

# CRUK control 8 samples
#pileupsC <- list.files("G:\\PolyA/faststorage/BACKUP/CRUK/plasma/N289", recursive = T, full.names = T, pattern = "bait.pileup")
pileupsC <- list.files("~/genomedk/PolyA/faststorage/BACKUP/CRUK/plasma/N289", recursive = T, full.names = T, pattern = "bait.pileup")
resCRUK <- lapply(pileupsC[1:8], FUN = function(x) sw_piles(pileup = pileupsC[1:8], pon = pon_obj2$pon, regions = pon_obj2$regions, prior = 0.5, model = "AND")) #Write the parameters instead of ...
names(resCRUK) <- basename(pileupsC[1:8])
saveRDS (resCRUK, "~/genomedk/matovanalysis/umiseq_analysis/CRUKcontrols/2020_11_06_CRUK_8CTL_SW.RDS")








adenoma_list <- info[ grepl("adenom", info$sample_type), "i" ] # 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 91
mahaA <- vector()
j=1
for (i in adenoma_list){
  mahaA[j] <- sum(countsW[i,,]* mafsC[i,,]/v1,  na.rm=T) 
    j = j + 1
}
plot(mahaA) # adenomas

control_list <- c(1,90, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104)
#<- length(unique(info[ grepl("control", info$sample_type), "i" ]))
mahaH <- vector()
j=1
for (i in control_list){
  mahaH[j] <- sum(countsW[i,,]* mafsC[i,,]/v1,  na.rm=T)
  j = j + 1
}
plot(mahaH) # controls

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


countsW <- readRDS("~/genomedk/matovanalysis/umiseq_analysis/R/cruki.RDS") # 

