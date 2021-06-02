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


# w3 per position
#w3pos <- mean (meanS - frli) * (meanS - frli)


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
########################
prior0 <- readRDS("~/genomedk/matovanalysis/umiseq_analysis/R/cosmic-prior.RDS") # 

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
  for (k in control_list[6:15]){
    mahaT[j] <- sum(countsW[i,,]* mafsC[k,,]/v1,  na.rm=T)/sum(countsW[i,,])#*prior0*prior0*prior0 #/sum(mafsC[k,,],  na.rm=T)
    print(sum(mafsC[k,,],  na.rm=T))
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
    mahaT2[j] <- sum(countsW[i,,]* mafsC[k,,]/v1,  na.rm=T)/sum(countsW[i,,])#*prior0*prior0*prior0 #/sum(mafsC[k,,],  na.rm=T)
    j = j + 1
  }
  #print(mahaT2)
}
length(mahaT2)
plot(mahaT2)

maha <- vector()
j=1
for (i in cancer_list){
  maha[j] <- sum(countsW[i,,]* mafsC[i,,]/v1,  na.rm=T) / sum(countsW[i,,])#*prior0*prior0*prior0 #/sum(mafsC[i,,],  na.rm=T)#45 of 95 cancers have WES, of them 37 have VAF>0
  #print(sum(countsW[i,,]))
  #print(sum(mafsC[i,,],  na.rm=T))
# maha[i] <- sum(countsW[6,,]* mafsC[i,,]/v1,  na.rm=T) # should be VAF mafsC[i,,]) #
if (maha[j]>0) {
j = j + 1
  }
#aux <- sum(countsW[6,,]*countsC001[i,,]/v1,  na.rm=T) # should be VAF mafsC[i,,]) #
# for countsC001 #6 [1] 142
# devide by PON var each pos
} 
print(maha)

mahaCancer <- maha[1:40] # 37 of 45 , w 8 cancers w VAFs 0
plot(mahaCancer) # Cancer samples with VAF>0

mahaControl <- c(mahaT)#, mahaT2) # 1575 (45x15+45x20)
plot(mahaT2,ylim=range(c(0,0.0035)))
plot(mahaT,ylim=range(c(0,0.0035)))
#ROC; find where overall success rate numbers are
rocTF = rep(0, (40+45*10))#15))#+45*20))
rocTF[1:40]=1

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

sort(mahaT,decreasing = T)[1:200]
#[1] 20.6824402 17.8332952 10.1667574  8.9166476  7.6151223  5.9444317  5.8567859  5.3905675  4.9414590  4.7634094  4.6975313  4.6848197
#[13]  4.6818134  4.4662639  3.9381400  3.9135748  3.7139431  3.2780695  3.2764601  3.2174697  3.1602292  3.1242025  3.0765785  2.9536050
#[25]  2.8877982  2.8804387  2.8606199  2.8506046  2.8346156  2.7498641  2.6952837  2.6952837  2.5935227  2.5873460  2.5685243  2.5643539
#[37]  2.5446014  2.3911765  2.3487657  2.2659081  2.2106057  2.2021961  2.1824921  2.1748940  2.1680396  2.1224018  2.1003417  2.0510523
#[49]  2.0209730  2.0205040  1.9280486  1.9115412  1.9078654  1.8581164  1.8358227  1.8350291  1.7906539  1.7893579  1.7728801  1.7377595
#[61]  1.7265585  1.7123495  1.7123495  1.7123495  1.6824171  1.6760842  1.5983281  1.5983281  1.5578680  1.5559973  1.5549756  1.5439182
#[73]  1.5278945  1.4708134  1.4677390  1.4639952  1.4611877  1.4454103  1.3311878  1.2921038  1.2853658  1.2842621  1.2658754  1.2469011
#[85]  1.2437829  1.2273125  1.2226587  1.2020833  1.2020833  1.2020833  1.2002153  1.1741323  1.1732490  1.1678085  1.1486386  1.1463694
#[97]  1.1447192  1.1400930  1.1279878  1.1198677  1.0879639  1.0714396  1.0649503  1.0587273  1.0221337  1.0104880  1.0104865  1.0032947
#[109]  1.0032947  0.9996725  0.9951214  0.9630789  0.9575009  0.9526133  0.9351531  0.9268816  0.9158949  0.9130055  0.9055762  0.8969626
#[121]  0.8947262  0.8829230  0.8826753  0.8797308  0.8602922  0.8602922  0.8569986  0.8514459  0.8442956  0.8422880  0.8404308  0.8064722
#[133]  0.7884524  0.7836316  0.7736735  0.7687767  0.7559319  0.7559319  0.7559319  0.7522368  0.7506946  0.7503197  0.7452887  0.7417887
#[145]  0.7354067  0.7290920  0.7267473  0.7250598  0.6981414  0.6880015  0.6851403  0.6810083  0.6768826  0.6709458  0.6640660  0.6625671
#[157]  0.6625671  0.6590524  0.6468919  0.6281224  0.6204087  0.6199989  0.6142378  0.6113293  0.5924450  0.5888188  0.5848648  0.5810856
#[169]  0.5669489  0.5639939  0.5504668  0.5499586  0.5439820  0.5266754  0.5191724  0.5167316  0.5020057  0.5011452  0.4997944  0.4969254
#[181]  0.4890635  0.4847036  0.4813434  0.4714797  0.4653496  0.4606719  0.4567602  0.4567602  0.4565028  0.4536893  0.4536893  0.4536893
#[193]  0.4511951  0.4435276  0.4363934  0.4311889  0.4311889  0.4275607  0.4179594  0.4179594

sort(mahaCancer)
#[1]  0.3666103  0.4513090  0.5001443  0.5110298  0.5161669  0.7061671  0.8449192  1.0352195  1.2850088  1.9004407  2.1350731  2.1500690
#[13]  2.2909422  2.9025490  3.3326619  3.3335790  3.8942093  4.9047943  5.3638151  5.8821141  5.9705000  6.2066127  6.4498850  7.2874429
#[25]  7.3366357  7.3555631  7.7998069  9.4741280  9.9544213 10.0973532 10.1381976 14.3441438 15.1550994 15.9179259 17.0618449 22.1524821
#[37] 60.8835776

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
hist(c(mahaT,mahaCancer),breaks = 60)
##############################SW##################################################

# CRUK control 8 samples
#pileupsC <- list.files("G:\\PolyA/faststorage/BACKUP/CRUK/plasma/N289", recursive = T, full.names = T, pattern = "bait.pileup")
pileupsC <- list.files("~/genomedk/PolyA/faststorage/BACKUP/CRUK/plasma/N289", recursive = T, full.names = T, pattern = "bait.pileup")
resCRUK <- lapply(pileupsC[1:8], FUN = function(x) sw_piles(pileup = pileupsC[1:8], pon = pon_obj2$pon, regions = pon_obj2$regions, prior = 0.5, model = "AND")) #Write the parameters instead of ...
names(resCRUK) <- basename(pileupsC[1:8])
saveRDS (resCRUK, "~/genomedk/matovanalysis/umiseq_analysis/CRUKcontrols/2020_11_06_CRUK_8CTL_SW.RDS")


aux = c( 5971 , 5310  ,4817 ,11835 , 6452 , 6409 , 8478 ,17439 , 9782 , 5547 ,16303,  7537 , 8103 , 8821,  8166,  7542,  9305 , 7983,  8699,  8551,  7443 ,17254 , 4431 , 4307 ,7843 ,2978 , 5093,  7637 , 3299,  3804 , 1154 , 3349 , 8559,  7092  ,2072,   887 , 9535 , 5891,  3527 , 3747,  9516  ,2896,  1489 , 2808 , 3124)





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

