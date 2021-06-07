library(stringr)
library(devtools)
library (LaplacesDemon)
library(FNN)
library("RColorBrewer")
library(ROCR)
library("ggplot2")
library("mixtools")
library(gplots)
library(spatstat)
library("transport")
library(readxl)
library(caret)
library(seewave)
library("glmnet")
library("ggpubr")
setwd ('G:\\DELFI_data/Derived/fragment_length_in_bins')
setwd ('~/genomedk/DELFI_data/Derived/fragment_length_in_bins')
require("reticulate")
py_install("pandas")
py_install("https://github.com/Hogfeldt/ctDNAtool")
source_python("~/genomedk/matovanalysis/DELFI_analysis/python/pickle_reader.py")
#pickle_data <- read_pickle_file("~/genomedk/DELFI2/Workspaces/per_and_elias/delfi2_length_5Mbp/DL000978HLQ0_100AM.pickle")
##excel_sheets(path = "~/genomedk/DELFI2/RawData/201217_Delfi2_fastq_and_sample_manifest_updated_batchinfo.xlsx")

#del2 <- read_excel("~/genomedk/DELFI2/RawData/201217_Delfi2_fastq_and_sample_manifest_updated_batchinfo.xlsx", sheet = 1) 
#end2 <- read_excel("~/genomedk/DELFI2/RawData/201217_Delfi2_fastq_and_sample_manifest_updated_batchinfo.xlsx", sheet = 2)# 1069 samples: 169 colon, 100 rectum, 800 control.
del2 <- read_excel("~/genomedk/DELFI2/Workspaces/matov/201217_Delfi2_fastq_and_sample_manifest_updated_batchinfo2.xlsx", sheet = 1) 
end2 <- read_excel("~/genomedk/DELFI2/Workspaces/matov/201217_Delfi2_fastq_and_sample_manifest_updated_batchinfo2.xlsx", sheet = 2)# 1069 samples: 169 colon, 100 rectum, 800 control.

col_list <- which(end2$diagnostic_group=="Colon cancer") # 169 samples
col1 <- which(end2$diagnostic_group=="Colon cancer" & end2$crc_stage=="I") # 19
col2 <- which(end2$diagnostic_group=="Colon cancer" & end2$crc_stage=="II") # 65
col3 <- which(end2$diagnostic_group=="Colon cancer" & end2$crc_stage=="III") # 42
col4 <- which(end2$diagnostic_group=="Colon cancer" & end2$crc_stage=="IV") # 43
col0 <- which(end2$diagnostic_group=="Adenoma colon" & end2$adenoma_risk=="HIGH") # 46 col Adenoma HIGH
colA <- which(end2$diagnostic_group=="Adenoma colon" & end2$adenoma_risk=="LOW") # 53 col Adenoma LOW
ctl1 <- which(end2$diagnostic_group=="No comorbidity-no finding") # 284
ctl2 <- which(end2$diagnostic_group=="Comorbidity-no finding") # 191
ctl3 <- which(end2$diagnostic_group=="Other finding") # 325
rec_list <- which(end2$diagnostic_group=="Rectal cancer") # 100 samples
rec1 <- which(end2$diagnostic_group=="Rectal cancer" & end2$crc_stage=="I") # 29
rec2 <- which(end2$diagnostic_group=="Rectal cancer" & end2$crc_stage=="II") # 25
rec3 <- which(end2$diagnostic_group=="Rectal cancer" & end2$crc_stage=="III") # 31
rec4 <- which(end2$diagnostic_group=="Rectal cancer" & end2$crc_stage=="IV") # 15
rec0 <- which(end2$diagnostic_group=="Adenoma rectum" & end2$adenoma_risk=="HIGH") # 21 rec Adenoma HIGH
recA <- which(end2$diagnostic_group=="Adenoma rectum" & end2$adenoma_risk=="LOW") # 14 rec Adenoma LOW

e2 <- data.frame(end2) # 1803   27 - stage and disease info comes from this file (sheet 2)
d2 <- data.frame(del2) # 765  10  - DELFI names (sheet 1)
m2 <- merge(e2,d2,by="SampleID") # 681  36 - only 681 are members of both groups

col_list2 <- which(m2$diagnostic_group=="Colon cancer") #  79
col12 <- which(m2$diagnostic_group=="Colon cancer" & m2$crc_stage=="I") #  8
col22 <- which(m2$diagnostic_group=="Colon cancer" & m2$crc_stage=="II") #  30
col32 <- which(m2$diagnostic_group=="Colon cancer" & m2$crc_stage=="III") #  18
col42 <- which(m2$diagnostic_group=="Colon cancer" & m2$crc_stage=="IV") #  23

rec12 <- which(m2$diagnostic_group=="Rectal cancer" & m2$crc_stage=="I") # 10
rec22 <- which(m2$diagnostic_group=="Rectal cancer" & m2$crc_stage=="II") # 11
rec32 <- which(m2$diagnostic_group=="Rectal cancer" & m2$crc_stage=="III") # 18
rec42 <- which(m2$diagnostic_group=="Rectal cancer" & m2$crc_stage=="IV") # 11

# 345 already analyzed, 7 colon replicates, 60 ctl3 later, 66 adenomas, tot. 478 for analysis (471 unique)
# plus 210 (of 681 in the list), all unique?, are "validation", i.e. m2$diagnostic_group is NA.
col02 <- which(m2$diagnostic_group=="Adenoma colon" & m2$adenoma_risk=="HIGH") # 20 col Adenoma HIGH
colA2 <- which(m2$diagnostic_group=="Adenoma colon" & m2$adenoma_risk=="LOW") # 28 col Adenoma LOW
rec02 <- which(m2$diagnostic_group=="Adenoma rectum" & m2$adenoma_risk=="HIGH") # 11 rec Adenoma HIGH
recA2 <- which(m2$diagnostic_group=="Adenoma rectum" & m2$adenoma_risk=="LOW") # 7 rec Adenoma LOW

sum(m2$Replicate) # 29
length(m2$Replicate)# 681
# these seem to be 681 Delfi IDs. 
val210 <- m2$diagnostic_group=="NA"
val210[is.na(val210)] <- 1 # list of 210 validation samples
val <- which(val210==1)
auxVAL <- unlist(sapply(m2$DELFI.ID[val] , function(x) grep(x, x = pileupsD2 )))
listVAL <- unique(auxVAL)

# 210 individual KLD to the 43 CTL Delfi1. pick top 189 bp . maybe go to 12 bins.
val189D2 <- valD2[,,189] #210 (samples) 574 (bins)

write.csv(t(val189D2[,1:555]),'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_val210_555bins_fr189.csv')
# 555 x 210 

dim(hnm)#  700 23310
hnm189 <- hnm[189,]#extract from 700 only 189 bp
h189x210 = array(0, dim=c(23310,210))
for (i in 1:210  ) {
h189x210[,i]<-hnm189 #duplicate the same vector 210 times to prepare the divergence.
}
write.csv(h189x210,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi1_healthy43_555bins_fr189.csv')
# 23310(42*555) x 210 
k189_val_d1ctl43 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2val210_D1ctl43_fr189.csv')
k189_vald1ctl43  <- k189_val_d1ctl43[2:211,2]
plot(k189_vald1ctl43 )
k189_val210_sorted<-sort(k189_vald1ctl43, index.return=TRUE,decreasing = T )
plot(k189_val210_sorted$x)
pileupsD2[listVAL[k189_val210_sorted$ix[1:21]]] # the colon cancer portion of the validation cohort
# $ix 26 184  78 200 188 149 148 116 103 115  74 157   9  85 109 134  47 112 182  16 205 
# plot vector of 210 divergence numbers in bits (its always for 189 bp only)
# are there peaks and dips?

# library no longer exists
classifier( method = c("randomForest", "svm", "nnet" ),  featureMat, positiveSamples, negativeSamples,  tunecontrol = tune.control(sampling = "cross", cross = 5), ...)


##for random forest, and using five-fold cross validation for obtaining optimal parameters
cl <- classifier( method = "randomForest", 
                      featureMat = hcTop20, 
                      positiveSamples = cv21, # col79x12 for frl195bp
                      negativeSamples = hv21,# ctl1_74x12 for frl195bp
                      tunecontrol = tune.control(sampling = "cross", cross = 5),
                      ntree = 100 ) #build 100 trees for the forest

# Error in family$linkfun(mustart) : Argument mu must be a nonempty numeric vector
simple_logistic_model = glm(data = data.frame(as.factor(df)), #dim(df) 153  12
                            family = binomial())
summary(simple_logistic_model)
# x 195 for COL, y 195 for CTL
glmnet(t(cv21), t(hv21), family = "binomial", alpha = 1, lambda = NULL)
#ROC; find where overall success rate numbers are
condition <- rbind(array(1, dim=c(79,12)), array(0, dim=c(74,12)))
pred <- prediction(hcTop20, condition, label.ordering = c(0, 1))  

#pearson 0.7627064 
cor(kc[1:499], kd2_col79, method = c("pearson", "kendall", "spearman"))
cor.test(kc[1:499], kd2_col79, method=c("pearson", "kendall", "spearman"))

plot(kc[1:499], ylim=range(c(0,2.7)), col="green", type="l")
par(new = TRUE)
plot(kd2_col79, ylim=range(c(0,2.7)), col="red", type="l")
legend(400,2.6,legend=c("DELFI1", "DELFI2"), col=c("green","red"),lty=1:1, cex=1.0)

plot(k_d2_col1, ylim=range(c(0,1.9)), col="green", type="l")
par(new = TRUE)
plot(k_d2_col2, ylim=range(c(0,1.9)), col="blue", type="l")
par(new = TRUE)
plot(k_d2_col3, ylim=range(c(0,1.9)), col="orange", type="l")
par(new = TRUE)
plot(k_d2_col4, ylim=range(c(0,1.9)), col="red", type="l")
legend(400,1.5,legend=c("CC SI", "CC SII", "CC SIII", "CC SIV"),col=c("green","blue","orange","red"),lty=1:1, cex=1.0)

plot(kd2_col79, ylim=range(c(0,.8)), col="green", type="l")
par(new = TRUE)
plot(kd2_col79_ctl2, ylim=range(c(0,.8)), col="blue", type="l")
par(new = TRUE)
plot(kd2_col79_ctl3, ylim=range(c(0,.8)), col="red", type="l")
legend(220,.8,legend=c("CC CTL1", "CC CTL2", "CC CTL3"),col=c("green","blue","red"),lty=1:1, cex=1.0)










m2$DELFI.ID[col_list2] # 79
#[1] "DL001860CRP0"   "DL001940CRP0"   "DL002182CRP0"   "DL002181CRP0"   "DL002023CRP0"   "DL001364CRP0_1" "DL001364CRP0"   "DL001503CRP0"   "DL001288CRP0"  
#[10] "DL001509CRP0"   "DL001422CRP0"   "DL001207CRP0"   "DL001800CRP0"   "DL001093CRP0"   "DL001785CRP0"   "DL001923CRP0"   "DL002241CRP0"   "DL001902CRP0"  
#[19] "DL001845CRP0"   "DL002194CRP0_1" "DL001578CRP0"   "DL001613CRP0"   "DL001614CRP0"   "DL001348CRP0"   "DL001096CRP0_1" "DL001096CRP0"   "DL001475CRP0"  
#[28] "DL001857CRP0"   "DL001835CRP0"   "DL001561CRP0"   "DL001142CRP0"   "DL001582CRP0_1" "DL001582CRP0"   "DL001448CRP0"   "DL002259CRP0"   "DL001581CRP0"  
#[37] "DL001888CRP0"   "DL002085CRP0"   "DL002000CRP0"   "DL002072CRP0"   "DL002267CRP0"   "DL001496CRP0"   "DL001390CRP0"   "DL001758CRP0_1" "DL001563CRP0"  
#[46] "DL002286CRP0"   "DL002286CRP0_1" "DL001849CRP0"   "DL002204CRP0"   "DL001550CRP0"   "DL001162CRP0"   "DL002227CRP0"   "DL001244CRP0"   "DL002280CRP0"  
#[55] "DL001236CRP0"   "DL002104CRP0"   "DL002046CRP0"   "DL002223CRP0"   "DL002123CRP0"   "DL001914CRP0"   "DL001918CRP0"   "DL001855CRP0"   "DL002026CRP0"  
#[64] "DL001270CRP0"   "DL001316CRP0"   "DL001562CRP0"   "DL001178CRP0"   "DL001572CRP0"   "DL001843CRP0"   "DL001843CRP0_1" "DL001964CRP0"   "DL001964CRP0_1"
#[73] "DL001763CRP0"   "DL001650CRP0"   "DL001965CRP0"   "DL002190CRP0_1" "DL002190CRP0"   "DL001793CRP0"   "DL001480CRP0"  

pileupsD2 <- list.files("~/genomedk/DELFI2/Workspaces/per_and_elias/delfi2_length_5Mbp", recursive = T, full.names = T, pattern = "tsv")
pileupsD2_1M <- list.files("~/genomedk/DELFI2/Workspaces/per_and_elias/delfi2_length_1Mbp", recursive = T, full.names = T, pattern = "tsv")

d2_test <- read.table(pileupsD2[204], header = TRUE)
d2_t2 <- read.table(pileupsD2[304], header = TRUE)
samplesREC <- sapply(m2$DELFI.ID[rec_list2] , function(x) grep(x, x = pileupsD2 )) 
samplesCOL <- sapply(m2$DELFI.ID[col_list2] , function(x) grep(x, x = pileupsD2 ))
listCOL1 <- sapply(m2$DELFI.ID[col12] , function(x) grep(x, x = pileupsD2 ))
listCOL2 <- sapply(m2$DELFI.ID[col22] , function(x) grep(x, x = pileupsD2 ))
listCOL3 <- unlist(sapply(m2$DELFI.ID[col32] , function(x) grep(x, x = pileupsD2 )))
listCOL4 <- unlist(sapply(m2$DELFI.ID[col42] , function(x) grep(x, x = pileupsD2 )))

listCOL0 <- unlist(sapply(m2$DELFI.ID[col02] , function(x) grep(x, x = pileupsD2 )))
listCOLA <- unlist(sapply(m2$DELFI.ID[colA2] , function(x) grep(x, x = pileupsD2 )))
listREC0 <- unlist(sapply(m2$DELFI.ID[rec02] , function(x) grep(x, x = pileupsD2 )))
listRECA <- unlist(sapply(m2$DELFI.ID[recA2] , function(x) grep(x, x = pileupsD2 )))
#######################################colon adenoma low##########################
col0D2 = array(0, dim=c(20,574,499))
j=1
for (i in 1:20 ) {
  print(i)
  #i=2
  auxFR <- read.table(pileupsD2[listCOL0[i]], header = TRUE) # sample per sample, file per file. 
  col0D2[j,,] <- unlist(auxFR[,2:500])
  j=j+1
}
colAD2 = array(0, dim=c(28,574,499))
j=1
for (i in 1:28 ) {
  print(i)
  #i=2
  auxFR <- read.table(pileupsD2[listCOLA[i]], header = TRUE) # sample per sample, file per file. 
  colAD2[j,,] <- unlist(auxFR[,2:500])
  j=j+1
}
rec0D2 = array(0, dim=c(11,574,499))
j=1
for (i in 1:11 ) {
  print(i)
  #i=2
  auxFR <- read.table(pileupsD2[listREC0[i]], header = TRUE) # sample per sample, file per file. 
  rec0D2[j,,] <- unlist(auxFR[,2:500])
  j=j+1
}
recAD2 = array(0, dim=c(7,574,499))
j=1
for (i in 1:7 ) {
  print(i)
  #i=2
  auxFR <- read.table(pileupsD2[listRECA[i]], header = TRUE) # sample per sample, file per file. 
  recAD2[j,,] <- unlist(auxFR[,2:500])
  j=j+1
}
col0D22 = array(0, dim=c(20*574,499))
colAD22 = array(0, dim=c(28*574,499))
rec0D22 = array(0, dim=c(11*574,499))
recAD22 = array(0, dim=c(7*574,499))
for (i in 1:499) {
  #i=1
  auxCOL0 <- col0D2[,,i]
  col0D22[,i] <- auxCOL0 
  auxCOLA<- colAD2[,,i]
  colAD22[,i] <- auxCOLA 
  #auxREC0 <- rec0D2[,,i]
  #rec0D22[,i] <- auxREC0 
  #auxRECA <- recAD2[,,i]
  #recAD22[,i] <- auxRECA 
}

write.csv(col0D22,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_col_adeH20.csv')
write.csv(colAD22,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_col_adeL28.csv')
write.csv(rec0D22,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_rec_adeH11.csv')
write.csv(recAD22,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_rec_adeL7.csv')

kd2_col0 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2col_adeH20_ctl1.csv')
k_d2_col0 <- kd2_col0[2:500,2]
plot(k_d2_col0)
kd2_colA <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2col_adeL28_ctl1.csv')
k_d2_colA <- kd2_colA[2:500,2]
plot(k_d2_colA)
kd2_rec0 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2rec_adeH11_ctl1.csv')
k_d2_rec0 <- kd2_rec0[2:500,2]
plot(k_d2_rec0)
kd2_recA <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2rec_adeL7_ctl1.csv')
k_d2_recA <- kd2_recA[2:500,2]
plot(k_d2_recA)
# combined figure for the four adenomas
plot(k_d2_col0, ylim=range(c(0,.36)) ,col="red", type="l")
par(new = TRUE)
plot(k_d2_colA, ylim=range(c(0,.36)), col="green", type="l")
legend(300,.3,legend=c("CC AH", "CC AL"),col=c("red","green"),lty=1:1, cex=1.0)

plot(k_d2_rec0, ylim=range(c(0,.7)), col="orange", type="l")
par(new = TRUE)
plot(k_d2_recA, ylim=range(c(0,.7)), col="blue", type="l")
legend(300,.6,legend=c( "RC AH", "RC AL"),col=c("orange","blue"),lty=1:1, cex=1.0)



auxCOL <- unlist(samplesCOL)

# auxCRC
#DL001860CRP0   DL001940CRP0   DL002182CRP0   DL002181CRP0   DL002023CRP0 DL001364CRP0_1  DL001364CRP01  DL001364CRP02   DL001503CRP0 
#377            419            511            510            452            188            188            189            250 
#DL001288CRP0   DL001509CRP0   DL001422CRP0   DL001207CRP0   DL001800CRP0   DL001093CRP0   DL001785CRP0   DL001923CRP0   DL002241CRP0 
#159            253            219            133            353             86            342            412            538 
#DL001902CRP0   DL001845CRP0 DL002194CRP0_1   DL001578CRP0   DL001613CRP0   DL001614CRP0   DL001348CRP0 DL001096CRP0_1  DL001096CRP01 
#398            370            519            277            292            293            184             87             87 
#DL001096CRP02   DL001475CRP0   DL001857CRP0   DL001835CRP0   DL001561CRP0   DL001142CRP0 DL001582CRP0_1  DL001582CRP01  DL001582CRP02 
#88            237            376            364            269            105            279            279            280 
#DL001448CRP0   DL002259CRP0   DL001581CRP0   DL001888CRP0   DL002085CRP0   DL002000CRP0   DL002072CRP0   DL002267CRP0   DL001496CRP0 
#226            542            278            391            471            442            467            544            248 
#DL001390CRP0 DL001758CRP0_1   DL001563CRP0  DL002286CRP01  DL002286CRP02 DL002286CRP0_1   DL001849CRP0   DL002204CRP0   DL001550CRP0 
#203            330            271            552            553            552            373            524            264 
#DL001162CRP0   DL002227CRP0   DL001244CRP0   DL002280CRP0   DL001236CRP0   DL002104CRP0   DL002046CRP0   DL002223CRP0   DL002123CRP0 
#110            531            143            549            139            480            460            529            490 
#DL001914CRP0   DL001918CRP0   DL001855CRP0   DL002026CRP0   DL001270CRP0   DL001316CRP0   DL001562CRP0   DL001178CRP0   DL001572CRP0 
#405            407            375            454            151            170            270            116            274 
#DL001843CRP01  DL001843CRP02 DL001843CRP0_1  DL001964CRP01  DL001964CRP02 DL001964CRP0_1   DL001763CRP0   DL001650CRP0   DL001965CRP0 
#368            369            368            429            430            429            333            305            431 
#DL002190CRP0_1  DL002190CRP01  DL002190CRP02   DL001793CRP0   DL001480CRP0 
#516            516            517            348            239 

listCOL <- unique(auxCOL)
#[1] 377 419 511 510 452 188 189 250 159 253 219 133 353  86 342 412 538 398 370 519 277 292 293 184  87  88 237 376 364 269 105 279 280 226 542 278 391 471 442 467 544 248 203 330
#[45] 271 552 553 373 524 264 110 531 143 549 139 480 460 529 490 405 407 375 454 151 170 270 116 274 368 369 429 430 333 305 431 516 517 348 239

length(unique(auxVAL)) # 210 (4 were replicated, even if its the same sample)
#################validation#######################
valD2 = array(0, dim=c(210,574,499))
j=1
for (i in 1:210  ) {
  print(i)
  #i=2
  auxFR <- read.table(pileupsD2[listVAL[i]], header = TRUE) # sample per sample, file per file. 
  valD2[j,,] <- unlist(auxFR[,2:500])
  j=j+1
}

length(unique(auxCOL)) # 79 (samples with indexes 188, 87, 279, 552, 368, 429, 516 were replicated, even if its the same sample)
#################colon all stages#######################
colD2 = array(0, dim=c(79,574,499))

#ctl1D2_1M = array(0, dim=c(74,2873,499)) # 5x574 genomic bins, i.e. 1Mb rather than 5Mb
j=1
for (i in 1:74 ) {#79
  print(i)
  #i=3
  auxFR <- read.table(pileupsD2[listCOL[i]], header = TRUE) # sample per sample, file per file. 
  #auxFR <- read.table(pileupsD2_1M[listCTL1[i]], header = TRUE) # sample per sample, file per file. 
  colD2[j,,] <- unlist(auxFR[,2:500])
  #ctl1D2_1M[j,,] <- unlist(auxFR[,2:500])
  
  # save 79 individual profiles as xls files
  #sa_name <- paste0('~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_col_individual', i,'.csv')
  print(length(unlist(auxFR[,2:500])))
  #auxSAVE <-matrix(unlist(auxFR[,2:500]),nrow=574,ncol=499)

  #write.csv(unlist(auxSAVE),sa_name)
  j=j+1
}
#######################################colon stage 1##########################
k_col1_i<-matrix(,nrow=8,ncol=700)
col1D2 = array(0, dim=c(8,574,499))
j=1
for (i in 1:8 ) {
  print(i)
  #i=3
  auxFR <- read.table(pileupsD2[listCOL1[i]], header = TRUE) # sample per sample, file per file. 
  col1D2[j,,] <- unlist(auxFR[,2:500])
  
  #sa_name <- paste0('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2col1_individual', i,'.csv')
  #aux<- read.csv(sa_name)
  #k_col1_i[i,] <- aux[2:500,2]
  j=j+1
}
#######################################colon stage 2##########################
col2D2 = array(0, dim=c(30,574,499))
j=1
for (i in 1:30 ) {
  print(i)
  #i=2
  auxFR <- read.table(pileupsD2[listCOL2[i]], header = TRUE) # sample per sample, file per file. 
  col2D2[j,,] <- unlist(auxFR[,2:500])
  j=j+1
}
#######################################colon stage 3##########################
col3D2 = array(0, dim=c(18,574,499))
j=1
for (i in 1:18 ) {
  print(i)
  #i=2
  auxFR <- read.table(pileupsD2[listCOL3[i]], header = TRUE) # sample per sample, file per file. 
  col3D2[j,,] <- unlist(auxFR[,2:500])
  j=j+1
}
#######################################colon stage 4##########################
col4D2 = array(0, dim=c(23,574,499))
j=1
for (i in 1:23 ) {
  print(i)
  #i=2
  auxFR <- read.table(pileupsD2[listCOL4[i]], header = TRUE) # sample per sample, file per file. 
  col4D2[j,,] <- unlist(auxFR[,2:500])
  j=j+1
}
##########################CONTROL##########################################################
ctl1_list <- which(m2$diagnostic_group=="No comorbidity-no finding") # 74

samplesCTL1 <- sapply(m2$DELFI.ID[ctl1_list] , function(x) grep(x, x = pileupsD2 )) 
auxCTL1 <- unlist(samplesCTL1)
#DL001791HLP0 DL001790HLP0 DL001437HLP0 DL001636HLP0 DL001481HLP0 DL001168HLP0 DL001171HLP0 DL001181HLP0 DL001215HLP0 DL001223HLP0 
#347          346          220          298          240          113          114          118          135          137 
#DL002195HLP0 DL002016HLP0 DL001609HLP0 DL001287HLP0 DL001440HLP0 DL001530HLP0 DL001092HLP0 DL001394HLP0 DL001839HLP0 DL002172HLP0 
#520          451          290          158          221          257           85          206          367          507 
#DL001775HLP0 DL001908HLP0 DL002199HLP0 DL002180HLP0 DL001605HLP0 DL001391HLP0 DL001884HLP0 DL001870HLP0 DL001872HLP0 DL001780HLP0 
#338          403          521          509          289          204          389          381          383          339 
#DL001956HLP0 DL001644HLP0 DL002088HLP0 DL001646HLP0 DL001798HLP0 DL001795HLP0 DL002031HLP0 DL001783HLP0 DL002074HLP0 DL001331HLP0 
#426          303          472          304          352          349          455          340          468          174 
#DL001451HLP0 DL001345HLP0 DL001382HLP0 DL001979HLP0 DL001553HLP0 DL001826HLP0 DL001868HLP0 DL001707HLP0 DL002129HLP0 DL002153HLP0 
#228          182          199          437          265          363          380          318          493          501 
#DL001148HLP0 DL002215HLP0 DL002053HLP0 DL001312HLP0 DL002050HLP0 DL001243HLP0 DL002206HLP0 DL001310HLP0 DL001272HLP0 DL002290HLP0 
#106          527          463          169          461          142          525          168          153          555 
#DL002008HLP0 DL002099HLP0 DL001408HLP0 DL001767HLP0 DL001279HLP0 DL002097HLP0 DL002092HLP0 DL001269HLP0 DL001457HLP0 DL001542HLP0 
#446          478          212          336          155          476          474          150          233          262 
#DL001961HLP0 DL001118HLP0 DL001117HLP0 DL001253HLP0 
#428           98           97          144 
listCTL1 <- unique(auxCTL1)

length(unique(auxCTL1)) # 74

ctl1D2 = array(0, dim=c(74,574,499))
#ctl1D2_1M = array(0, dim=c(74,2873,499)) # 5x574 genomic bins, i.e. 1Mb rather than 5Mb
j=1
for (i in 1:74  ) {
  print(i)
  #i=2
  auxFR <- read.table(pileupsD2[listCTL1[i]], header = TRUE) # sample per sample, file per file. 
  
  ctl1D2[j,,] <- unlist(auxFR[,2:500])
  
  # save 74 individual profiles as xls files
  #sa_name <- paste0('~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_ctl1_individual', i,'.csv')
  #auxSAVE <-matrix(unlist(auxFR[,2:500]),nrow=574,ncol=499)
  
  
  #write.csv(unlist(auxSAVE),sa_name)
  j=j+1
}

ctl1D22 = array(0, dim=c(74*574,499))
col1D22 = array(0, dim=c(8*574,499))
col2D22 = array(0, dim=c(30*574,499))
col3D22 = array(0, dim=c(18*574,499))
col4D22 = array(0, dim=c(23*574,499))

colD22 = array(0, dim=c(79*574,499))

ctl1D22_1M = array(0, dim=c(74*2873,499))
colD22_1M = array(0, dim=c(79*2873,499))
for (i in 1:499) {
  #i=1
  auxCTL <- ctl1D2[,,i]
  ctl1D22[,i] <- auxCTL 
  auxCOL <- colD2[,,i]
  colD22[,i] <- auxCOL
  #auxCOL1 <- col1D2[,,i]
  #col1D22[,i] <- auxCOL1 
  #auxCOL2<- col2D2[,,i]
  #col2D22[,i] <- auxCOL2 
  #auxCOL3 <- col3D2[,,i]
  #col3D22[,i] <- auxCOL3 
  auxCOL4 <- col4D2[,,i]
  col4D22[,i] <- auxCOL4 
}
write.csv(colD22_1M,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_col_all79_1M.csv')
write.csv(ctl1D22_1M,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_ctl1_74_1M.csv')

k_d2_col79_1M <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2_COL79_CTL1_1M.csv')
kd2_col79_1M <- k_d2_col79_1M[2:575,2]
plot(kd2_col79_1M)

k_d2_crc129 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2_CRC129_ctl1.csv')
kd2_crc129 <- k_d2_crc129[2:575,2]
plot(kd2_crc129)

ctlD22 <- rbind(ctl1D22,ctl2D22)
ctlD22 <- rbind(ctlD22,ctl3D22)
write.csv(ctlD22,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_ctl_all276.csv')
k_d2_crc129_ctl276 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2_CRC129_ctl276.csv')
kd2_crc129_ctl276 <- k_d2_crc129_ctl276[2:575,2]
plot(kd2_crc129_ctl276)
k_d2_col79_ctl276 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2_COL79_ctl276.csv')
kd2_col79_ctl276 <- k_d2_col79_ctl276[2:575,2]
plot(kd2_col79_ctl276)




write.csv(colD22,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_col_all79.csv')

dim(colD2)# 79 574 499
dim(ctl1D2)# 74 574 499
c195<-colD2[,,195]
h195<-ctl1D2[,,195]
hh195<-ctl11D2[,,195]
dim(c195) # 79 574
dim(h195)#  74 574
dim(hh195)#  73 574
write.csv(c195,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_col_fr195.csv')
write.csv(h195,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_ctl1_fr195.csv')
write.csv(hh195,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_ctl73_fr195.csv')

kk195_d2_col <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2_COLfr195_ctl73.csv')
kk195d2_col <- kk195_d2_col[2:575,2]
plot(kk195d2_col)
indx<-which(kk195d2_col>-5)#.847)#  141 153 158 170 218 347 374 405 408 458 528 558
cv22 <-matrix(,nrow=79,ncol=length(indx))
cv22<-c195[,indx]
hv22 <-matrix(,nrow=73,ncol=length(indx))
hv22<-hh195[,indx]
hcTop20 <- rbind(cv22 , hv22)
df<-scale(hcTop20)
df[is.nan(df)] <- 0
col <- colorRampPalette(brewer.pal(11, "RdYlBu"))(256)
hm <- heatmap(df, scale = "none", col =  col) 
k195_d2_col <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2_COLfr195_ctl1.csv')
k195d2_col <- k195_d2_col[2:575,2]
plot(k195d2_col)

ind195<-which(k195d2_col>-5)#  all 
ind195<-which(k195d2_col>1.12)#  33  48  90 145 156 191 213 356 429 483 529 
ind195<-which(k195d2_col>1.72)#   356  
ind195<-which(k195d2_col>1.58)#   356  529
indx<-ind195
cv21 <-matrix(,nrow=79,ncol=length(indx))
cv21<-c195[,indx]
hv21 <-matrix(,nrow=74,ncol=length(indx))
hv21<-h195[,indx]

hcTop20 <- rbind(cv21 , hv21)
df<-scale(hcTop20)
col <- colorRampPalette(brewer.pal(11, "RdYlBu"))(256)
hm <- heatmap(df, scale = "none", col =  col) 


c195<-colD2[,,195]
h195<-ctl1D2[,,195]

cv2N <- matrix(colD22[,195],nrow=79,ncol=574)
hv2N <- matrix(ctl1D22[,195],nrow=74,ncol=574)
hcTop20 <- rbind(cv2N , hv2N)#(c195,h195)#(cv2N , hv2N)
df<-scale(hcTop20)
col <- colorRampPalette(brewer.pal(11, "RdYlBu"))(256)
hm <- heatmap(df, scale = "none", col =  col) 






c365<-colD2[,,365]
h365<-ctl1D2[,,365]
write.csv(c365,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_col_fr365.csv')
write.csv(h365,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_ctl1_fr365.csv')

k365_d2_col <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2_COLfr365_ctl1.csv')
k365d2_col <- k365_d2_col[2:575,2]
plot(k365d2_col)

indx<-which(k365d2_col>-5)#.847)#  141 153 158 170 218 347 374 405 408 458 528 558
cv22 <-matrix(,nrow=79,ncol=length(indx))
cv22<-c365[,indx]
hv22 <-matrix(,nrow=74,ncol=length(indx))
hv22<-h365[,indx]

hcTop20 <- rbind(cv22 , hv22)
df<-scale(hcTop20)
df[is.nan(df)] <- 0
col <- colorRampPalette(brewer.pal(11, "RdYlBu"))(256)
hm <- heatmap(df, scale = "none", col =  col) 

cv2<-t(rbind(t(cv21),t(cv22)))
hv2<-t(rbind(t(hv21),t(hv22)))

hcTop20 <- rbind(cv2 , hv2)
df<-scale(hcTop20)
col <- colorRampPalette(brewer.pal(11, "RdYlBu"))(256)
hm <- heatmap(df, scale = "none", col =  col) 

hnm499<-hnm[1:499,]
cnm499<-cnm[1:499,]
write.csv(t(cnm499),'~/genomedk/matovanalysis/DELFI_analysis/python/delfi1_crc27_frl499.csv')
write.csv(t(hnm499),'~/genomedk/matovanalysis/DELFI_analysis/python/delfi1_ctl43_frl499.csv')

k1_d2_col79 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2col79_D1ctl43.csv')
k1d2_col79 <- k1_d2_col79[2:500,2]
plot(k1d2_col79)

k1_d2_col79hg38 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2col79_D1ctl43HG38.csv')
k1d2_col79hg38 <- k1_d2_col79hg38[2:500,2]
plot(k1d2_col79hg38)

k2_d2_col79hg38 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2col79_D1ctl86HG38.csv')
k2d2_col79hg38 <- k2_d2_col79hg38[2:500,2]
plot(k2d2_col79hg38)

hM189 <- matrix(hnm499[189,], ncol = 555, byrow = 42) #  
write.csv(hM189,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi1_ctl43_frl499_189.csv')

cM189 <- matrix(colD22[,189], ncol = 574, byrow = 79) #  
write.csv(cM189[,1:555],'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_col79_bin555_189.csv')

col_bins = array(0, dim=c(79*499,574))
ctl1_bins = array(0, dim=c(74*499,574))
#col_bins <- matrix(, ncol = 574, byrow = 79) #  
dim(colD2) # 79 574 499
for (i in 1:574){
  auxCOL <- colD2[,i,]
  col_bins[,i] <- auxCOL
  auxCTL1 <- ctl1D2[,i,]
  ctl1_bins[,i] <- auxCTL1
}
#ctl1_bins <- matrix(ctl1D22, ncol = 574, byrow = 74) # 

write.csv(ctl1_bins,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_ctl1Bins.csv')
write.csv(col_bins,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_colBins.csv')
k_d2_colBins <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2col_Bins.csv')
kd2_colBins<-k_d2_colBins[2:701,2]
plot(kd2_colBins)


dim(colD22) # 45346 499
dim(ctl1D22) # 42476   499
colD2B10<- array(0, dim = c(79*574*10,50))
ctl1D2B10<- array(0, dim = c(74*574*10,50))
umicB10<-array(0,dim=c(67*595*10,50))
for (i in 0:49){
  #i=2
  k=10*i+1
  if (i != 49) {
    colD2B10[,i] <- colD22[,k:(k+9)]
    ctl1D2B10[,i] <- ctl1D22[,k:(k+9)]
  } else {
    colD2B10[,i] <- colD22[,(k-1):(k+8)]# repeat one FRl from the previous segment
    ctl1D2B10[,i] <- ctl1D22[,(k-1):(k+8)]# repeat one FRl from the previous segment
  }
}
write.csv(colD2B10,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_col_79_binned10frl.csv')
write.csv(ctl1D2B10,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_ctl1_74_binned10frl.csv')



k1_d2_col79_189 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2col79_D1ctl43_fr189.csv')
k1d2_col79_189 <- k1_d2_col79_189[2:556,2]
plot(k1d2_col79_189)

indx<-which(k1d2_col79_189>-5)#1.24) #1.3
cv22 <-matrix(,nrow=79,ncol=length(indx))
cv22<-cM189[,indx]
hv22 <-matrix(,nrow=42,ncol=length(indx))
hv22<-hM189[,indx]

hcTop20 <- rbind(cv22 , hv22)
df<-scale(hcTop20)
col <- colorRampPalette(brewer.pal(11, "RdYlBu"))(256)
hm <- heatmap(df, scale = "none", col =  col) 

hgA <- hist(colD22[,195], breaks = 50 , plot = FALSE) # Save first histogram data
hgB <- hist(ctl1D22[,195],breaks = 50, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,470), ylim = c(0,3500)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,470), ylim = c(0,3500)) # Add 2nd histogram using different color

hgA <- hist(colD22[,365], breaks = 30 , plot = FALSE) # Save first histogram data
hgB <- hist(ctl1D22[,365],breaks = 30, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,70), ylim = c(0,5000)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,70), ylim = c(0,5000)) # Add 2nd histogram using different color

hgA <- hist(colD22[,137], breaks = 150 , plot = FALSE) # Save first histogram data
hgB <- hist(ctl1D22[,137],breaks = 50, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,520), ylim = c(0,4200)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,520), ylim = c(0,4200)) # Add 2nd histogram using different color

hgA <- hist(col4D22[,137], breaks = 150 , plot = FALSE) # Save first histogram data
hgB <- hist(ctl1D22[,137],breaks = 50, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,520), ylim = c(0,4200)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,520), ylim = c(0,4200)) # Add 2nd histogram using different color

hgA <- hist(colD22[,277], breaks = 100 , plot = FALSE) # Save first histogram data
hgB <- hist(ctl1D22[,277],breaks = 30, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,40), ylim = c(0,6000)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,40), ylim = c(0,6000)) # Add 2nd histogram using different color

hgA <- hist(col4D22[,277], breaks = 100 , plot = FALSE) # Save first histogram data
hgB <- hist(ctl1D22[,277],breaks = 30, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,40), ylim = c(0,6000)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,40), ylim = c(0,6000)) # Add 2nd histogram using different color



plot(kc4[1:499], ylim=range(c(0,2.4)), col="red", type="l")
par(new = TRUE)
plot(kd2_col4, ylim=range(c(0,2.4)), col="green", type="l")
par(new = TRUE)
plot(k2d2_col4, ylim=range(c(0,2.4)), col="blue", type="l")
legend(240,2.48,legend=c("D1 8 CRC Stage IV 43 CTL", "D2 23 CC Stage IV 74 CTL No comorbidity", "D2 23 CC Stage IV 71 CTL Comorbidity"),col=c("red","green","blue"),lty=1:1, cex=1.0)

#pearson  
cor(kc4[1:499], kd2_col4, method = c("pearson", "kendall", "spearman"))
cor(kc4[1:499], k2d2_col4, method=c("pearson", "kendall", "spearman"))
cor(kd2_col4, k2d2_col4, method=c("pearson", "kendall", "spearman"))
max_k2d2_colStage <- c(max(k2d2_col1),max(k2d2_col2),max(k2d2_col3),max(k2d2_col4))
d2_maxSens_spec90 <- c(.62,.52,.52,.97)
cor(max_k2d2_colStage, d2_maxSens_spec90, method=c("pearson", "kendall", "spearman"))
d2_maxSens_spec95 <- c(.50,.42,.43,.96)
cor(max_kd2_colStage, d2_maxSens_spec95, method=c("pearson", "kendall", "spearman"))
cor(max_kd2_colStage, d2_maxSens_spec90, method=c("pearson", "kendall", "spearman"))
max_kd2_colStage <- c(max(kd2_col1),max(kd2_col2),max(kd2_col3),max(kd2_col4))

k_d2_col_i1 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2col_individual1.csv')
kd2_col_i[1,] <- k_d2_col_i1[2:500,2]
plot(kd2_col_i[1,])
k_d2_col_i2 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2col_individual2.csv')
kd2_col_i[2,] <- k_d2_col_i2[2:500,2]
plot(kd2_col_i[2,])

k_d2_col_i3 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2col_individual3.csv')
kd2_col_i[3,]<- k_d2_col_i3[2:500,2]
plot(kd2_col_i[3,])
k_d2_col_i4 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2col_individual4.csv')
kd2_col_i[4,] <- k_d2_col_i4[2:500,2]
plot(kd2_col_i[4,])

kd2_col_i<-matrix(,nrow=79,ncol=499)
for(i in 1:79) {
  sa_name <- paste0('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2col_individual', i,'.csv')
  aux <- read.csv(sa_name)
  kd2_col_i[i,] <- aux[2:500,2]
  print(which.max(aux[2:500,2]))
  #plot(kd2_col_i[i,], ylim=range(c(0,2.4)), col="red", type="l"))  
  #par(new = TRUE)
}
for(i in 1:79) {
  #i=4
#plot(kd2_col_i[i,], ylim=range(c(0,5)),col= rgb(1,0,0,i/10)) # looks good
  plot(kd2_col_i[i,], ylim=range(c(0,5)),col= col2rgb(i, alpha = FALSE)) # looks good
  
  par(new = TRUE)
}

col_list2 <- which(m2$diagnostic_group=="Colon cancer") #  79
colS1<-kd2_col_i[m2$crc_stage[col_list2]=="I",]
colS2<-kd2_col_i[m2$crc_stage[col_list2]=="II",]
colS3<-kd2_col_i[m2$crc_stage[col_list2]=="III",]
colS4<-kd2_col_i[m2$crc_stage[col_list2]=="IV",]
mS1=NULL
mS2=NULL
mS3=NULL
mS4=NULL
maS1=NULL
maS2=NULL
maS3=NULL
maS4=NULL
for(i in 1:8) {
  plot(colS1[i,], ylim=range(c(0,5.2)), col= col2rgb(i+50, alpha = FALSE), type="l")    
  mS1[i]<-which.max(colS1[i,])
  maS1[i]<-max(colS1[i,])
  par(new = TRUE)
}
#plot(mS1, maS1, pch = 24, cex=2, col="blue", bg="red", lwd=2)# plot the dot of the max
#par(new = TRUE)

for(i in 1:30) {
  plot(colS2[i,], ylim=range(c(0,5.2)), col= col2rgb(2*i+50), type="l") 
  mS2[i]<-which.max(colS2[i,])
  maS2[i]<-max(colS2[i,])
  par(new = TRUE)
}
for(i in 1:18) {
  plot(colS3[i,], ylim=range(c(0,5.2)), col= col2rgb(2*i+50), type="l")
  mS3[i]<-which.max(colS3[i,])
  maS3[i]<-max(colS3[i,])
  par(new = TRUE)
}
for(i in 1:23) {
  plot(colS4[i,], ylim=range(c(0,5.2)), col= col2rgb(2*i+50), type="l") 
  mS4[i]<-which.max(colS4[i,])
  maS4[i]<-max(colS4[i,])
  par(new = TRUE)
}

dim(colS1)# 8 499
dim(colS2)#30 499
dim(colS3)#   18 499
dim(k_crc4_i)#   23 499
#dim(k_ctl1_i)# 43 700
colS <- rbind(colS1,colS2)
colS <- rbind(colS,colS3)
colS <- rbind(colS,colS4)

c198<-colS[,70] 
c364<-colS[,253] 
h198<-k_ctl1_i[,70] 
h364<-k_ctl1_i[,253] 
plot(c198,c364,xlim=c(0,1.3),ylim=c(0,2),col="red")
par(new = TRUE)
plot(h198, h364,xlim=c(0,1.3),ylim=c(0,2),col="green")
frl1=207; frl2=366; k=5.4
# boundary for capacity regions for Stage I, II, III, IV.
plot(colS1[,frl1],colS1[,frl2],xlim=c(0,k),ylim=c(0,k),col="blue",pch=1)
par(new = TRUE)
plot(colS2[,frl1],colS2[,frl2],xlim=c(0,k),ylim=c(0,k),col="orange",pch=2)
par(new = TRUE)
plot(colS3[,frl1],colS3[,frl2],xlim=c(0,k),ylim=c(0,k),col="red",pch=3)
par(new = TRUE)
plot(colS4[,frl1],colS4[,frl2],xlim=c(0,k),ylim=c(0,k),col="brown",pch=4)
par(new = TRUE)
plot(kd2_ctl1_i[,frl1],kd2_ctl1_i[,frl2],xlim=c(0,k),ylim=c(0,k),col="green",pch=5)
legend(3.5,2,legend=c("KLD for Col SI samples", "KLD for Col SII samples", "KLD for Col SIII samples","KLD for Col SIV samples","KLD for Ctl1 74 samples"),col=c("blue","orange","red","brown","green"),lty=1:1, cex=1.0)

plot(abs(colSums(colS)-colSums(kd2_ctl1_i)))
which.max(abs(colSums(colS,na.rm = T)-colSums(kd2_ctl1_i,na.rm = T)))
which(abs(colSums(colS)-colSums(kd2_ctl1_i))>50)







hgA <- hist(mS1, breaks = 50 , plot = FALSE) # Save first histogram data
hgB <- hist(mS2,breaks = 50, plot = FALSE) # Save 2nd histogram data
hgC <- hist(mS3, breaks = 50 , plot = FALSE) # Save first histogram data
hgD <- hist(mS4,breaks = 50, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(100,400), ylim = c(0,8)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(100,400), ylim = c(0,8)) # Add 2nd histogram using different color
plot(hgC, col = rgb(1,0,1,1/10),add = TRUE,xlim = c(100,400), ylim = c(0,8)) # Plot 1st histogram using a transparent color
plot(hgD, col = rgb(0,1,1,1/10), add = TRUE,xlim = c(100,400), ylim = c(0,8)) # Add 2nd histogram using different color

plot(mS1, ylim=range(c(100,400)),col="green")
par(new = TRUE)
plot(mS2, ylim=range(c(100,400)), col="blue")
par(new = TRUE)
plot(mS3, ylim=range(c(100,400)), col="orange")
par(new = TRUE)
plot(mS4, ylim=range(c(100,400)), col="red")
legend(240,2.48,legend=c("D1 8 CRC Stage IV 43 CTL", "D2 23 CC Stage IV 74 CTL No comorbidity", "D2 23 CC Stage IV 71 CTL Comorbidity"),col=c("red","green","blue"),lty=1:1, cex=1.0)



k_d2_ctl1_i1 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2ctl1_individual1.csv')
kd2_ctl1_i[1,] <- k_d2_ctl1_i1[2:500,2]
plot(kd2_ctl1_i[1,])
k_d2_ctl1_i2 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2ctl1_individual2.csv')
kd2_ctl1_i[2,] <- k_d2_ctl1_i2[2:500,2]
plot(kd2_ctl1_i[2,])

# brew a palette with 74 values and pass it to "col" with an index i.
kd2_ctl1_i<-matrix(,nrow=74,ncol=499)
mC1=NULL
maC1=NULL
for(i in 1:74) {# break down the plot of 74 into 4 plots of 18
  #i=1
  sa_name <- paste0('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2ctl1_individual', i,'.csv')
  aux <- read.csv(sa_name)
  kd2_ctl1_i[i,] <- aux[2:500,2]
  t1<-floor(i/9)
  k <- (i-t1*9)
  plot(kd2_ctl1_i[i,],  ylim=range(c(0,3.8)), col= palette()[k], type="l") # there is something wrong with the colors
  print(palette()[k])#print(col2rgb(2*i+50))
  mC1[i]<-which.max(kd2_ctl1_i[i,] )
  maC1[i]<-max(kd2_ctl1_i[i,] )
  par(new = TRUE)
}
print(median(maC1))
print(median(mC1))

dim(k_ctl1_i) # 43 700 # 9 is NA
totKLD_D1_CTL1 <- rowSums(k_ctl1_i[,1:499], dims = 1)# take only 1:499
dim(kd2_ctl1_i) #74 499
totKLD_D2_CTL1 <- rowSums(kd2_ctl1_i, dims = 1)# is it much higher? 
hgA <- hist(totKLD_D1_CTL1, breaks = 10, plot = FALSE) # Save first histogram data
hgB <- hist(totKLD_D2_CTL1,  breaks = 30, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,1000), ylim = c(0,24)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,1000), ylim = c(0,24)) # Add 2nd histogram using different color

ctl1D2[35,,]
totFRs_D2_CTL1 <- rowSums(ctl1D2)
totFRs_D2_CTL1[35] # 38186686 outlier
mean(totFRs_D2_CTL1)+3*sd(totFRs_D2_CTL1) # 33584769
plot(totFRs_D1_CTL1)
plot(totFRs_D2_CTL1)
shapiro.test(totFRs_D1_CTL1)
shapiro.test(totFRs_D2_CTL1)
totFRs_D2_CTL11<-vector()
totFRs_D2_CTL11[1:34]<-totFRs_D2_CTL1[1:34]
totFRs_D2_CTL11[35:73]<-totFRs_D2_CTL1[36:74]
shapiro.test(totFRs_D2_CTL11)
totKLD_D1_CTL11 <- rowSums(k_ctl1_i, dims = 1)# take only 1:499
plot(totFRs_D1_CTL11)
shapiro.test(totFRs_D1_CTL11)

c_195<-colD2[,,195]
h1_195<-ctl1D2[,,195]
c2M195 <- matrix(c_195, ncol = 574, byrow = 79) #  
h2M195 <- matrix(h1_195, ncol = 574, byrow = 74) # 
mean(h2M195[,539])/mean(c2M195[,539])# GPX4  1.362526
mean(h2M195[,356])/mean(c2M195[,356]) # VDAC2 1.314116
hgA <- hist(c2M195[,356], breaks = 30 , plot = FALSE) # Save first histogram data
hgB <- hist(h2M195[,356],breaks = 20, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,400), ylim = c(0,10)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,400), ylim = c(0,10)) # Add 2nd histogram using different color

c_205<-colD2[,,205]
h1_205<-ctl1D2[,,205]
c2M205 <- matrix(c_205, ncol = 574, byrow = 79) #  
h2M205 <- matrix(h1_205, ncol = 574, byrow = 74) # 
mean(h2M205[,539])/mean(c2M205[,539])# GPX4 1.394561
hgA <- hist(c2M205[,363], breaks = 30 , plot = FALSE) # Save first histogram data
hgB <- hist(h2M205[,363],breaks = 20, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,270), ylim = c(0,13)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,270), ylim = c(0,13)) # Add 2nd histogram using different color

# Stage I 205 bp, Stage II 137 bp, Stage III 364 bp, Stage IV 364 bp   
c1_205<-col1D2[,,205]
h1_205<-ctl1D2[,,205]
write.csv(c1_205,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_col1_fr205.csv')
write.csv(h1_205,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_ctl1_fr205.csv')
#######################
col277<-colD2[,,277]
h1_277<-ctl1D2[,,277]
cM277 <- matrix(col277, ncol = 574, byrow = 79) #  
hM277 <- matrix(h1_277, ncol = 574, byrow = 74) # 
hgA <- hist(cM364[,406], breaks = 20 , plot = FALSE) # Save first histogram data
hgB <- hist(hM364[,406],breaks = 40, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,60), ylim = c(0,5)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,60), ylim = c(0,5)) # Add 2nd histogram using different color
##
hgA <- hist(cM277[,426], breaks = 50 , plot = FALSE) # Save first histogram data
hgB <- hist(hM277[,426],breaks = 20, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,25), ylim = c(0,15)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,25), ylim = c(0,15)) # Add 2nd histogram using different color
#######################
col365<-colD2[,,365]
h1_365<-ctl1D2[,,365]
cM365 <- matrix(col365, ncol = 574, byrow = 79) #  
hM365 <- matrix(h1_365, ncol = 574, byrow = 74) # 
##
hgA <- hist(cM365[,405], breaks = 50 , plot = FALSE) # Save first histogram data
hgB <- hist(hM365[,405],breaks = 50, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,60), ylim = c(0,10)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,60), ylim = c(0,10)) # Add 2nd histogram using different color

c4_137<-col4D2[,,137]
write.csv(c4_137,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_col4_fr137.csv')
k137_d2_col4 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2col4_fr137_ctl1.csv')
k137d2_col4 <- k137_d2_col4[2:575,2]
plot(k137d2_col4)

h1_277<-ctl1D2[,,277]
write.csv(h1_277,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_ctl1_fr277.csv')
c4_277<-col4D2[,,277]
cM277 <- matrix(c4_277, ncol = 574, byrow = 23) #  
hM277 <- matrix(h1_277, ncol = 574, byrow = 74) # 
hgA <- hist(cM277[,426], breaks = 50 , plot = FALSE) # Save first histogram data
hgB <- hist(hM277[,426],breaks = 5, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,88), ylim = c(0,22)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,88), ylim = c(0,22)) # Add 2nd histogram using different color
mean(cM277[,426])/mean(hM277[,426]) # 

write.csv(c4_277,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_col4_fr277.csv')
k277_d2_col4 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2col4_fr277_ctl1.csv')
k277d2_col4 <- k277_d2_col4[2:575,2]
plot(k277d2_col4)
#######################
c2_137<-col2D2[,,137]
h1_137<-ctl1D2[,,137]
cM137 <- matrix(c4_137, ncol = 574, byrow = 23) #  
hM137 <- matrix(h1_137, ncol = 574, byrow = 74) # 
hgA <- hist(cM137[,555], breaks = 50 , plot = FALSE) # Save first histogram data
hgB <- hist(hM137[,555],breaks = 25, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,400), ylim = c(0,10)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,400), ylim = c(0,10)) # Add 2nd histogram using different color
mean(cM137[,555])/mean(hM137[,555]) #1.785717

write.csv(c2_137,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_col2_fr137.csv')
write.csv(h1_137,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_ctl1_fr137.csv')
k137_d2_col2 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2col2_fr137_ctl1.csv')
k137d2_col2 <- k137_d2_col2[2:575,2]
plot(k137d2_col2)
indx<-which(k137d2_col2>0.33)#12
cv <-matrix(,nrow=30,ncol=length(indx))
cv<-c2_137[,indx]
hv <-matrix(,nrow=74,ncol=length(indx))
hv<-h1_137[,indx]
hcTop20 <- rbind(cv , hv)
df<-scale(hcTop20)
df[is.nan(df)] <- 0
col <- colorRampPalette(brewer.pal(11, "RdYlBu"))(256)
hm <- heatmap(df, scale = "none", col =  col) 

c3_364<-col3D2[,,364]
h1_364<-ctl1D2[,,364]
write.csv(c3_364,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_col3_fr364.csv')
write.csv(h1_364,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_ctl1_fr364.csv')
k364_d2_col3 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2col3_fr364_ctl1.csv')
k364d2_col3 <- k364_d2_col3[2:575,2]
plot(k364d2_col3)
indx<-which(k364d2_col3>0.08)#12
cv <-matrix(,nrow=18,ncol=length(indx))
cv<-c3_364[,indx]
hv <-matrix(,nrow=74,ncol=length(indx))
hv<-h1_364[,indx]
hcTop20 <- rbind(cv , hv)
df<-scale(hcTop20)
df[is.nan(df)] <- 0
col <- colorRampPalette(brewer.pal(11, "RdYlBu"))(256)
hm <- heatmap(df, scale = "none", col =  col) 

c4_364<-col3D2[,,364]
write.csv(c4_364,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_col4_fr364.csv')
k364_d2_col4 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2col4_fr364_ctl1.csv')
k364d2_col4 <- k364_d2_col4[2:575,2]
plot(k364d2_col4)
indx<-which(k364d2_col4>0.08)#12
cv <-matrix(,nrow=23,ncol=length(indx))
cv<-c4_364[,indx]
hv <-matrix(,nrow=74,ncol=length(indx))
hv<-h1_364[,indx]
hcTop20 <- rbind(cv , hv)
df<-scale(hcTop20)
df[is.nan(df)] <- 0
col <- colorRampPalette(brewer.pal(11, "RdYlBu"))(256)
hm <- heatmap(df, scale = "none", col =  col) 
########################################################################################
write.csv(col1D22,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_col1_all8.csv')
write.csv(col2D22,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_col2_all30.csv')
write.csv(col3D22,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_col3_all18.csv')
write.csv(col4D22,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_col4_all23.csv')

write.csv(ctl1D22,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_ctl1_74.csv')
ctl11D22= array(0, dim=c(73*574,499))
ctl11D2= array(0, dim=c(73,574,499))
ctl11D2[1:34,,]<-ctl1D2[1:34,,]
ctl11D2[35:73,,]<-ctl1D2[36:74,,]
for (i in 1:499) {
  #i=1
  auxCTL11 <- ctl11D2[,,i]
  ctl11D22[,i] <- auxCTL11 
}
write.csv(ctl11D22,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_ctl1_73.csv')
k2_d2_col79 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2_COL79_ctl73.csv')
k2d2_col79 <- k2_d2_col79[2:500,2]
plot(k2d2_col79)
cor(kd2_col79, k2d2_col79, method=c("pearson", "kendall", "spearman"))
cor(kc[1:499], k2d2_col79, method=c("pearson", "kendall", "spearman"))
#cor(kc[1:499], kd2_col79, method=c("pearson", "kendall", "spearman"))

k_d2_col79 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2_COL79.csv')
kd2_col79 <- k_d2_col79[2:500,2]
plot(kd2_col79)
##############COLON STAGE 1 vs CTL1
k_d2_col1 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2_COL1_8_ctl1.csv')
kd2_col1 <- k_d2_col1[2:500,2]
plot(kd2_col1)
##############COLON STAGE 2 vs CTL1
k_d2_col2 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2_COL2_30_ctl1.csv')
kd2_col2 <- k_d2_col2[2:500,2]
plot(kd2_col2)
##############COLON STAGE 3 vs CTL1
k_d2_col3 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2_COL3_18_ctl1.csv')
kd2_col3 <- k_d2_col3[2:500,2]
plot(kd2_col3)
##############COLON STAGE 4 vs CTL1
k_d2_col4 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2_COL4_23_ctl1.csv')
kd2_col4 <- k_d2_col4[2:500,2]
plot(kd2_col4)
##############COLON STAGE 1 vs CTL2
k2_d2_col1 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2_COL1_8_ctl2.csv')
k2d2_col1 <- k2_d2_col1[2:500,2]
plot(k2d2_col1)
##############COLON STAGE 2 vs CTL2
k2_d2_col2 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2_COL2_30_ctl2.csv')
k2d2_col2 <- k2_d2_col2[2:500,2]
plot(k2d2_col2)
##############COLON STAGE 3 vs CTL2
k2_d2_col3 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2_COL3_18_ctl2.csv')
k2d2_col3 <- k2_d2_col3[2:500,2]
plot(k2d2_col3)
##############COLON STAGE 4 vs CTL2
k2_d2_col4 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2_COL4_23_ctl2.csv')
k2d2_col4 <- k2_d2_col4[2:500,2]
plot(k2d2_col4)


##########################COLON REVERSED KLD############################
k2_d2_col79 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2_COL79_revered.csv')
k2d2_col79 <- k2_d2_col79[2:500,2]
plot(k2d2_col79)
#########################################################################################
ctl2_list <- which(m2$diagnostic_group=="Comorbidity-no finding") # 71

samplesCTL2 <- sapply(m2$DELFI.ID[ctl2_list] , function(x) grep(x, x = pileupsD2 )) 
auxCTL2 <- unlist(samplesCTL2)

listCTL2 <- unique(auxCTL2)

length(unique(auxCTL2)) # 71

ctl2D2 = array(0, dim=c(71,574,499))
j=1
for (i in 1:71  ) {
  print(i)
  #i=2
  auxFR <- read.table(pileupsD2[listCTL2[i]], header = TRUE) # sample per sample, file per file. 
  ctl2D2[j,,] <- unlist(auxFR[,2:500])
  j=j+1
}
ctl2D22 = array(0, dim=c(71*574,499))

for (i in 1:499) {
  #i=1
  auxCTL2 <- ctl2D2[,,i]
  ctl2D22[,i] <- auxCTL2 
}
write.csv(ctl2D22,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_ctl2_71.csv')

k_d2_col79_ctl2 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2_COL79_ctl2_71.csv')
kd2_col79_ctl2 <- k_d2_col79_ctl2[2:500,2]
plot(kd2_col79_ctl2)
##################################################################################
ctl3_list <- which(m2$diagnostic_group=="Other finding") # 131

samplesCTL3 <- sapply(m2$DELFI.ID[ctl3_list] , function(x) grep(x, x = pileupsD2 )) 
auxCTL3 <- unlist(samplesCTL3)

listCTL3 <- unique(auxCTL3)

length(unique(auxCTL3)) # 131

ctl3D2 = array(0, dim=c(131,574,499)) # super computer runs out of memory for 131 samples. 
j=1
for (i in 1:131  ) {
  print(i)
  #i=2
  auxFR <- read.table(pileupsD2[listCTL3[i]], header = TRUE) # sample per sample, file per file. 
  ctl3D2[j,,] <- unlist(auxFR[,2:500])
  j=j+1
}
ctl3D22 = array(0, dim=c(131*574,499)) # super computer runs out of memory for 131 samples. 

for (i in 1:499) {
  #i=1
  auxCTL3 <- ctl3D2[,,i]
  ctl3D22[,i] <- auxCTL3 
}
write.csv(ctl3D22,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_ctl3_71.csv') # the first 71 of 131

k_d2_col79_ctl3 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2_COL79_ctl3_71.csv')
kd2_col79_ctl3 <- k_d2_col79_ctl3[2:500,2]
plot(kd2_col79_ctl3)
##################################
rec_list2 <- which(m2$diagnostic_group=="Rectal cancer") # 

m2$DELFI.ID[rec_list2] # 50

samplesREC <- sapply(m2$DELFI.ID[rec_list2] , function(x) grep(x, x = pileupsD2 )) 
auxREC <- unlist(samplesREC)

listREC <- unique(auxREC)

length(unique(auxREC)) # 79 (samples with indexes 188, 87, 279, 552, 368, 429, 516 were replicated, even if its the same sample)

recD2 = array(0, dim=c(50,574,499))
j=1
for (i in 1:50  ) {
  print(i)
  #i=2
  auxFR <- read.table(pileupsD2[listREC[i]], header = TRUE) # sample per sample, file per file. 
  recD2[j,,] <- unlist(auxFR[,2:500])
  j=j+1
}
recD22 = array(0, dim=c(50*574,499))
for (i in 1:499) {
  #i=1
  auxREC <- recD2[,,i]
  recD22[,i] <- auxREC 
}
write.csv(recD22,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_rec_all50.csv')
k_d2_rec50 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2_REC50_ctl1.csv')
kd2_rec50 <- k_d2_rec50[2:500,2]
plot(kd2_rec50)

write.csv(rbind(recD22,colD22),'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_crc_all129.csv')
k_d2_rec50 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2_REC50_ctl1.csv')
kd2_rec50 <- k_d2_rec50[2:500,2]
plot(kd2_rec50)


#################################################################################################
k2_d2_rec50 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2_REC50_ctl2.csv')
k2d2_rec50 <- k2_d2_rec50[2:500,2]
plot(k2d2_rec50)
#################################################################################################
k3_d2_rec50 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2_REC50_ctl3.csv')
k3d2_rec50 <- k3_d2_rec50[2:500,2]
plot(k3d2_rec50)
#################################################################################################

plot(fft(kd2_col79, inverse = FALSE))
kccF <- fft(kd2_col79, inverse = FALSE)
eigen(kccF)# expects a matrix as an input; add as many kld vectors as the number of samples.
# compute kld spectrum followed by sample classification based on the spectrum
plot(fft(kc, inverse = FALSE))
plot(fft(kd2_rec50, inverse = FALSE))
plot(fft(k2d2_rec50, inverse = FALSE))
plot(fft(k3d2_rec50, inverse = FALSE))
plot(abs(fft(kc, inverse = FALSE)))
plot(abs(fft(kd2_col79, inverse = FALSE)))
plot(abs(fft(kd2_rec50, inverse = FALSE)))

# fft(kc, inverse = FALSE) should be taken square because of negative values; zero component should be removed




#bFr <- read.delim("filtered_window_data_5MB_1_700_m5000.txt")
#bf <- read.table("filtered_window_data_5MB_1_700_m5000.txt")
bFr <- read.table("~/genomedk/DELFI_data/Derived/fragment_length_in_bins/filtered_window_data_5MB_1_700_m5000.txt", header = TRUE)
typeof(bFr) # list
#lookUp <- read.table("5MB_map_file.txt")

# originial paper has 504 bins per patient sample, which covers 2.52G of genome
s1 = unlist(bFr[1,3:702]) # 700 (fragments) length histogram values for one 5mps bin
# up to 49 bins per chromosome, some get very few
s2 = unlist(bFr[2,3:702])
hist(s1[s1>100])
hist(s2[s2>100])

#sL <- read.csv("U:\\Documents/R/delfi_scripts-master/sample_reference.csv")
#sL <- read.csv("~/genomedk/matovanalysis/DELFI_analysis/R/sample_reference.csv")
sL <- read.csv("~/genomedk/DELFI_data/sample_reference.csv")

iL <- sL[sL[,6]=="Lung Cancer",3] # list  
iD <- sL[sL[,6]=="Cholangiocarcinoma",3] # list  
iB <- sL[sL[,6]=="Breast Cancer",3] # list  
iG <- sL[sL[,6]=="Gastric cancer",3] # list  
iC <- sL[sL[,6]=="Colorectal Cancer",3] # list of CRC patient numbers/names,
iC1 <- sL[sL[,6]=="Colorectal Cancer" & sL[,7]=="I",3] # list of CRC patient numbers/names,
iC2 <- sL[sL[,6]=="Colorectal Cancer" & sL[,7]=="II",3] # list of CRC patient numbers/names,
iC3 <- sL[sL[,6]=="Colorectal Cancer" & sL[,7]=="III",3] # list of CRC patient numbers/names,
iC4 <- sL[sL[,6]=="Colorectal Cancer" & sL[,7]=="IV",3] # list of CRC patient numbers/names,
iO <- sL[sL[,6]=="Ovarian Cancer",3] # list of ovarian cancer patient numbers/names,
iP <- sL[sL[,6]=="Pancreatic Cancer",3] # list  
iH <- sL[sL[,6]=="Healthy",3] # list of healthy/control samples, 215
iU  <- sL[sL[,6]=="Duodenal Cancer",3] # list  

length(bFr[bFr[,1]==iC[1],1])# 555 bins, i.e number of rows with CRC FRs for sample 1
# number of fragments per patient
sum(crc1[,3:702])#54026559

write.csv(cccD22,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_ccc195.csv')

k_d2_ccc195_ctl276 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2_CCC195_ctl276.csv')
kd2_ccc195_ctl276   <- k_d2_ccc195_ctl276  [2:500,2]
plot(kd2_ccc195_ctl276 )

names(colD22) <- names(recD22)
cc1D22 <- rbind (colD22, recD22)
names(recAD22) <- names(cc1D22)
cc2D22 <- rbind (cc1D22, recAD22)
names(rec0D22) <- names(cc2D22)
cc3D22 <- rbind (cc2D22, rec0D22)
names(colAD22) <- names(cc3D22)
cc4D22 <- rbind (cc3D22, colAD22)
names(col0D22) <- names(cc4D22)
cccD22 <- rbind (cc4D22, col0D22)


write.csv(ctl,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi1_HG38_ctl215.csv')
write.csv(ccc,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi1_HG38_ccc208.csv')

k_d1_ccc208_ctl215 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD1HG38_CCC208_ctl215.csv')
kd1_ccc208_ctl215   <- k_d1_ccc208_ctl215 [2:500,2]
plot(kd1_ccc208_ctl215 )

names(dcc) <- names(lcc)
cc1 <- rbind (dcc, lcc)
names(bcc) <- names(cc1)
cc2 <- rbind (cc1, bcc)
names(gcc) <- names(cc2)
cc3 <- rbind (cc2, gcc)
names(crc) <- names(cc3)
cc4 <- rbind (cc3, crc)
names(ovc) <- names(cc4)
cc5 <- rbind (cc4, ovc)
names(pcc) <- names(cc5)
cc6 <- rbind (cc5, pcc)
names(ucc) <- names(cc5)
ccc <- rbind (cc6, ucc)

listD1UCC1 <- unlist(sapply(iU, function(x) grep(x, x = pileupsD1 )))
ucc <- read.table(pileupsD1[listD1UCC1], header = TRUE)[,2:500]

listD1LCC12 <- unlist(sapply(iL, function(x) grep(x, x = pileupsD1 )))
nbLCC <- length(iL)
lcc <-matrix(,nrow=574,ncol=500)#nrow=555,ncol=702)
# for all patients
for(i in 1:nbLCC) {
  #i= 1
  #lcc1 <- bFr[bFr[,1]==iL[i],]
  print(i)
  auxFR <- read.table(pileupsD1[listD1LCC12[i]], header = TRUE) # sample per sample, file per file. 
  lcc1 <- auxFR[,2:500]
  print(sum(lcc1)) # 555 x 702   
  if (i==1){
    lcc <- lcc1 
  } else {
    names(lcc1) <- names(lcc)
    lcc <- rbind(lcc , lcc1) 
  }
}

listD1DCC25 <- unlist(sapply(iD, function(x) grep(x, x = pileupsD1 )))
nbDCC <- length(iD)
dcc <-matrix(,nrow=574,ncol=500)#nrow=555,ncol=702)
# for all patients
for(i in 1:nbDCC) {
  #i= 1
  #dcc1 <- bFr[bFr[,1]==iD[i],]
  print(i)
  auxFR <- read.table(pileupsD1[listD1DCC25[i]], header = TRUE) # sample per sample, file per file. 
  dcc1 <- auxFR[,2:500]
  print(sum(dcc1)) # 555 x 702 
  if (i==1){
    dcc <- dcc1 
  } else {
    names(dcc1) <- names(dcc)
    dcc <- rbind(dcc , dcc1) 
  }
}

listD1BCC54 <- unlist(sapply(iB, function(x) grep(x, x = pileupsD1 )))
nbBCC <- length(iB)
bcc <-matrix(,nrow=574,ncol=500)#nrow=555,ncol=702)
# for all patients
for(i in 1:nbBCC) {
  #i= 1
  #bcc1 <- bFr[bFr[,1]==iB[i],]
  print(i)
  auxFR <- read.table(pileupsD1[listD1BCC54[i]], header = TRUE) # sample per sample, file per file. 
  bcc1 <- auxFR[,2:500]
  print(sum(bcc1)) # 555 x 702     
  if (i==1){
    bcc <- bcc1 
  } else {
    names(bcc1) <- names(bcc)
    bcc <- rbind(bcc , bcc1) 
  }
}

listD1GCC27 <- unlist(sapply(iG, function(x) grep(x, x = pileupsD1 )))
nbGCC <- length(iG)
gcc <-matrix(,nrow=574,ncol=500)#nrow=555,ncol=702)
# for all patients
for(i in 1:nbGCC) {
  #i= 1
  #gcc1 <- bFr[bFr[,1]==iG[i],]
  print(i)
  auxFR <- read.table(pileupsD1[listD1GCC27[i]], header = TRUE) # sample per sample, file per file. 
  gcc1 <- auxFR[,2:500]
  print(sum(gcc1)) # 555 x 702    
  if (i==1){
    gcc <- gcc1 
  } else {
    names(gcc1) <- names(gcc)
    gcc <- rbind(gcc , gcc1) 
  }
}

listD1PCC34 <- unlist(sapply(iP, function(x) grep(x, x = pileupsD1 )))
pcc <-matrix(,nrow=574,ncol=500)#nrow=555,ncol=702)
nbPCC <- length(iP)
# for all patients
for(i in 1:nbPCC) {
  #i= 1
  #pcc1 <- bFr[bFr[,1]==iP[i],]
  print(i)
  auxFR <- read.table(pileupsD1[listD1PCC34[i]], header = TRUE) # sample per sample, file per file. 
  pcc1 <- auxFR[,2:500]
  print(sum(pcc1)) # 555 x 702  
  if (i==1){
    pcc <- pcc1 
  } else {
    names(pcc1) <- names(pcc)
    pcc <- rbind(pcc , pcc1)#  
  }
}

listD1OVC28 <- unlist(sapply(iO, function(x) grep(x, x = pileupsD1 )))
ovc <-matrix(,nrow=574,ncol=500)#nrow=555,ncol=702)
nbOVC <- length(iO)
# for all patients
for(i in 1:nbOVC) {
  #i= 1
  print(i)
  #ovc1 <- bFr[bFr[,1]==iO[i],]
  auxFR <- read.table(pileupsD1[listD1OVC28[i]], header = TRUE) # sample per sample, file per file. 
  ovc1 <- auxFR[,2:500]
  print(sum(ovc1)) # 555 x 702
  if (i==1){
    ovc <- ovc1 
  } else {
    names(ovc1) <- names(ovc)
    ovc <- rbind(ovc , ovc1) 
  }
}





# UMIIMPROVE
pileupsUI <- list.files("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/umiseq_paper/divergence/data/57PRE5Mb", recursive = T, full.names = T, pattern = "tsv")
nbUMII <- length(pileupsUI)
umiiD2 <- array(0,c(nbUMII,595,499))
umii1 <-matrix(,nrow=595,ncol=499)#nrow=555,ncol=702)
j=1
for(i in 1:nbUMII) {
  #i= 2
  print(i)
  auxFR <- read.table( pileupsUI[i], header = TRUE) # sample per sample, file per file. 
  testU <- as.integer(unlist(auxFR[,2:500]))
  umii1 <- t(matrix(testU, ncol = dim(auxFR)[1], byrow = (dim(auxFR)[2]-2)) )# convert back to matrix form
  umiiD2[j,,] <- umii1
  j=j+1
  print(dim(umii1)) # 555 x 702
  if (i==1){
    umii <- umii1
  } else {
    names(umii) <- names(umii1) 
    umii <- rbind(umii , umii1)
  }
}
plot(rowSums(umii1))
plot(colSums(umii1))

umiiM <- matrix(umii, ncol = 595, byrow = 499) #

# compare WASSERSTEIN for umiseq CRUK, IMPROVE and Delfi1 CRC and Delfi2 COL
wa_umii<-vector()
wa_umic<-vector()
for (i in 1:499){
  wa_umii[i] <- wasserstein1d(umii[,i],umi[,i])
  wa_umic[i] <- wasserstein1d(umic[,i],umi[,i])
}





umid <- rbind(umii, umic)

dim(umii) # 17850x499 for 30PreOps, 33320x499 for 56 PreOps
hist(testU2[364,], breaks = 150, ylim = c(0, 10))
hist(umii[,364], breaks = 150, ylim = c(0, 200))# dense up to 200 on x axis
hist(umii[,198], breaks = 150, ylim = c(0, 200))# dense up to 500 on x axis
# UMICRUK
pileupsUC <- list.files("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/umiseq_paper/divergence/data/CRUK5Mb", recursive = T, full.names = T, pattern = "tsv")
umic<-vector()
nbUMIC <- length(pileupsUC)
umicD2 <- array(0,c(nbUMIC,595,499))
umic1 <-matrix(,nrow=595,ncol=499)#nrow=555,ncol=702)
j=1
for(i in 1:nbUMIC) {
  #i= 1
  print(i)
  auxFR <- read.table( pileupsUC[i], header = TRUE) # sample per sample, file per file. 
  testU <- as.integer(unlist(auxFR[,2:500]))
  umic1 <- t(matrix(testU, ncol = dim(auxFR)[1], byrow = (dim(auxFR)[2]-2)) )# convert back to matrix form
  umicD2[j,,] <- umic1
  j=j+1
  print(dim(umic1)) # 555 x 702
  if (i==1){
    umic <- umic1
  } else {
    names(umic) <- names(umic1) 
    umic <- rbind(umic , umic1)
  }
}
dim(umic) # 40460 x 499 for  68 CRUK PreOps
#UMISEQ
auxFR <- read.table( "~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/umiseq_paper/divergence/data/length_matrix1.tsv", header = TRUE) # sample per sample, file per file. 
#testS <- as.integer(unlist(auxFR[,2:500]))
#testS2 <- matrix(testS, ncol = dim(auxFR)[1], byrow = (dim(auxFR)[2]-2)) # convert back to matrix form
pileupsU <- list.files("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/umiseq_paper/divergence/data/45PON5Mb", recursive = T, full.names = T, pattern = "tsv")
nbUMI <- length(pileupsU)
umiD2 <- array(0,c(nbUMI,595,499))
umi1 <-matrix(,nrow=574,ncol=499)#nrow=555,ncol=702)
j=1
for(i in 1:nbUMI) {
  #i= 2
  print(i)
  auxFR <- read.table( pileupsU[i], header = TRUE) # sample per sample, file per file. 
  testS <- as.integer(unlist(auxFR[,2:500]))
  umi1 <- t(matrix(testS, ncol = dim(auxFR)[1], byrow = (dim(auxFR)[2]-2)) )# convert back to matrix form
  umiD2[j,,] <- umi1
  j=j+1
  print(dim(umi1)) # 555 x 702
  if (i==1){
    umi <- umi1
  } else {
    names(umi) <- names(umi1) 
    umi <- rbind(umi , umi1)
  }
}

dim(umi) # 17850   499 for 30PONs, 26775   499 for 45 PONs
hist(testS2[364,], breaks = 150, ylim = c(0, 10))
hist(umi[,364], breaks = 150, ylim = c(0, 200))# dense up to 200 on x axis
hist(umi[,198], breaks = 150, ylim = c(0, 200),xlim = c(0,800))# dense up to 500 on x axis

hgA <- hist(umii[,198],breaks = 150 , plot = FALSE) # Save first histogram data
hgB <- hist(umi[,198] ,breaks = 150, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,800), ylim = c(0,150)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,800), ylim = c(0,150)) # Add 2nd histogram using different color


hgA <- hist(umii[,356],breaks = 250 , plot = FALSE) # Save first histogram data
hgB <- hist(umi[,356] ,breaks = 150, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,300), ylim = c(0,400)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,300), ylim = c(0,400)) # Add 2nd histogram using different color







hgA <- hist(cv130[,363],breaks = 50 , plot = FALSE) # Save first histogram data
hgB <- hist(hv130[,363] ,breaks = 50, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,55), ylim = c(0,4)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,55), ylim = c(0,4)) # Add 2nd histogram using different color


hgA <- hist(cv356[,505],breaks = 40 , plot = FALSE) # Save first histogram data
hgB <- hist(hv356[,505] ,breaks = 20, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,200), ylim = c(0,8)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,200), ylim = c(0,8)) # Add 2nd histogram using different color

write.csv(umi,'~/genomedk/matovanalysis/DELFI_analysis/python/umiseq_pon45.csv')
write.csv(umii,'~/genomedk/matovanalysis/DELFI_analysis/python/umiseq_pre56.csv')
write.csv(umic,'~/genomedk/matovanalysis/DELFI_analysis/python/umiseq_CRpre68.csv')
write.csv(umid,'~/genomedk/matovanalysis/DELFI_analysis/python/umiseq_totalPre124.csv')

k_umi_dPre124 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceUMIseq_totalPre124_Pon45.csv')
k_umidPre124  <- k_umi_dPre124[2:500,2]
plot(k_umidPre124)

k_umi_cPre68 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceUMIseqCRPre68Pon45.csv')
k_umicPre68  <- k_umi_cPre68[2:500,2]
plot(k_umicPre68)
# fit k_umicPre67
x <- seq(1, 499, by=1)
#fit4 <- lm(k_umicPre67~poly(x))
fit4 <- lm(k_umicPre67~poly(x,4,raw=TRUE))
xx <- seq(1,499, length=499)
plot(x,k_umicPre67,pch=19)
lines(xx, predict(fit4, data.frame(x=xx)), col="red")

fitTest <- rbind(k_umi_cPre67, x)
plot(k_umicPre67, cex = 0.5)
hpts <- chull(fitTest, y=NULL)
hpts <- c(hpts, hpts[1])
lines(k_umicPre67[hpts ])

X <- array(0, dim=c(2,499))
chull(X)
## Not run: 
# Example usage from graphics package
plot(X, cex = 0.5)
hpts <- chull(X)
hpts <- c(hpts, hpts[1])
lines(X[hpts, ])


umi199 <- umi[,199]
hv199 <- matrix(umi199, ncol = 595, byrow = 45) #

umic199 <- umic[,199]
cv199 <- matrix(umic199, ncol = 595, byrow = 68) #

write.csv(hv199,'~/genomedk/matovanalysis/DELFI_analysis/python/umiseq_pon45_frl199.csv')
write.csv(cv199,'~/genomedk/matovanalysis/DELFI_analysis/python/umiseq_CRpre68_frl199.csv')

k_umi_cPre68_frl199 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceUMIseqCRPre68Pon45_frl199.csv')
k_umicPre68_frl199  <- k_umi_cPre68_frl199[2:500,2]
plot(k_umicPre68_frl199)

indC <- which(k_umicPre68_frl199>0.2)

hcTop20 <- rbind(cv199[,indC], hv199[,indC])
df<-scale(hcTop20)
col <- colorRampPalette(brewer.pal(11, "RdYlBu"))(256)
hm <- heatmap(df, scale = "none", col =  col) 

k_umi_iPre30 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceUMIseqPre30Pon30.csv')
k_umiiPre30  <- k_umi_iPre30[2:500,2]
plot(k_umiiPre30)

k_umi_iPre56 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceUMIseqPre56Pon45.csv')
k_umiiPre56  <- k_umi_iPre56[2:500,2]
plot(k_umiiPre56)

paramEsts= gmdistribution.fit(k_umiiPre56,499)  

mixmdl = normalmixEM(k_umiiPre56)
plot(mixmdl,which=2)
lines(density(k_umiiPre56), lty=2, lwd=2)

tab <- data.frame(x=k_umiiPre56, r=499)
#Apply function nls
res <- nls( r ~ k*exp(-1/2*(x-mu)^2/sigma^2), start=c(mu=15,sigma=5,k=1) , data = tab)




x <- seq(1, 499, by=1)
fit4 <- lm(k_umiiPre56~poly(x,4,raw=TRUE))
xx <- seq(1,499, length=499)
plot(x,k_umiiPre56,pch=19)
lines(xx, predict(fit4, data.frame(x=xx)), col="red")


x <- seq(1, 499, by=1)
fitTest <- t(rbind(x,k_umicPre68))
plot(fitTest)
hpts <- chull(fitTest)
hpts <- c(hpts, hpts[1])
lines(fitTest[hpts, ])
# CONVEX HULL
x <- seq(1, 499, by=1)
fitTest <- t(rbind(x,k_umiiPre56))
plot(fitTest)
hpts <- chull(fitTest)
hpts <- c(hpts, hpts[1])
lines(fitTest[hpts, ])

plot(k_umiiPre56)



umi130 <- umi[,130]
hv130 <- matrix(umi130, ncol = 595, byrow = 30) #

umii130 <- umii[,130]
cv130 <- matrix(umii130, ncol = 595, byrow = 30) #

write.csv(hv130,'~/genomedk/matovanalysis/DELFI_analysis/python/umiseq_pon30_frl130.csv')
write.csv(cv130,'~/genomedk/matovanalysis/DELFI_analysis/python/umiseq_pre30_frl130.csv')

k_umi_iPre30_f130 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceUMIseqPre30Pon30_frl130.csv')
k_umiiPre30_f130  <- k_umi_iPre30_f130[2:596,2]
plot(k_umiiPre30_f130)
indU <- which(k_umiiPre30_f130>0.366)
# with most values
umi199 <- umi[,199]
hv199 <- matrix(umi199, ncol = 595, byrow = 45) #
umic199 <- umic[,199]
cv199 <- matrix(umic199, ncol = 595, byrow = 68) #CRUK

######################################
CRUK68 = array(0, dim=c(68,595,499))
PON45 = array(0, dim=c(45,595,499))
for (i in 1:499 ) {
  print(i)
  umic_aux <- umic[,i]
  CRUK68[,,i] <- matrix(umic_aux, ncol = 595, byrow = 68) #CRUK
  
  umi_aux <- umi[,i]
  PON45[,,i] <- matrix(umi_aux, ncol = 595, byrow = 45) # PON
}

#####################
umi361 <- umi[,361]
hv361 <- matrix(umi361, ncol = 595, byrow = 45) #
umic361 <- umic[,361]
cv361 <- matrix(umic361, ncol = 595, byrow = 68) #CRUK

hgA <- hist(cv361[,220],breaks = 6 , plot = FALSE) # Save first histogram data
hgB <- hist(hv361[,220] ,breaks = 6, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,130), ylim = c(0,7)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,130), ylim = c(0,7)) # Add 2nd histogram using different color



hgA <- hist(cv199[,303],breaks = 1 , plot = FALSE) # Save first histogram data
hgB <- hist(hv199[,303] ,breaks = 2, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,2), ylim = c(0,70)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,2), ylim = c(0,70)) # Add 2nd histogram using different color

hgA <- hist(cv199[,341],breaks = 2 , plot = FALSE) # Save first histogram data
hgB <- hist(hv199[,341] ,breaks = 4, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,3), ylim = c(0,60)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,3), ylim = c(0,60)) # Add 2nd histogram using different color

hgA <- hist(cv199[,60],breaks = 1 , plot = FALSE) # Save first histogram data
hgB <- hist(hv199[,60] ,breaks = 2, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,2), ylim = c(0,66)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,2), ylim = c(0,66)) # Add 2nd histogram using different color


umi356 <- umi[,356]
hv356 <- matrix(umi356, ncol = 595, byrow = 45) #
umii356 <- umii[,356]
cv356 <- matrix(umii356, ncol = 595, byrow = 56) #

hgA <- hist(cv356[,209],breaks = 50 , plot = FALSE) # Save first histogram data
hgB <- hist(hv356[,209] ,breaks = 30, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,110), ylim = c(0,7)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,110), ylim = c(0,7)) # Add 2nd histogram using different color

hgA <- hist(cv356[,220],breaks = 50 , plot = FALSE) # Save first histogram data
hgB <- hist(hv356[,220] ,breaks = 30, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,150), ylim = c(0,7)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,150), ylim = c(0,7)) # Add 2nd histogram using different color

write.csv(hv356,'~/genomedk/matovanalysis/DELFI_analysis/python/umiseq_pon45_frl356.csv')
write.csv(cv356,'~/genomedk/matovanalysis/DELFI_analysis/python/umiseq_pre56_frl356.csv')

k_umi_iPre56_f356 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceUMIseqPre56Pon45_frl356.csv')
k_umiiPre56_f356  <- k_umi_iPre56_f356[2:596,2]
plot(k_umiiPre56_f356)
indU2 <- which(k_umiiPre56_f356>0.5)

hcTop20 <- rbind(cv356[,indU2], hv356[,indU2])
df<-scale(hcTop20)
col <- colorRampPalette(brewer.pal(11, "RdYlBu"))(256)
hm <- heatmap(df, scale = "none", col =  col) 


hgA <- hist(cv130[,196],breaks = 50 , plot = FALSE) # Save first histogram data
hgB <- hist(hv130[,196] ,breaks = 20, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,4000), ylim = c(0,4)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,4000), ylim = c(0,4)) # Add 2nd histogram using different color

hgA <- hist(cv130[,363],breaks = 50 , plot = FALSE) # Save first histogram data
hgB <- hist(hv130[,363] ,breaks = 50, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,55), ylim = c(0,4)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,55), ylim = c(0,4)) # Add 2nd histogram using different color

hcTop20 <- rbind(cv130[,indU], hv130[,indU])
df<-scale(hcTop20)
col <- colorRampPalette(brewer.pal(11, "RdYlBu"))(256)
hm <- heatmap(df, scale = "none", col =  col) 


dim(k_crc1_i)# 7 700
dim(k_crc2_i)#7 700
dim(k_crc3_i)#   5 700
dim(k_crc4_i)#   8 700
dim(k_ctl1_i)# 43 700
k_crc_i <- rbind(k_crc1_i,k_crc2_i)
k_crc_i <- rbind(k_crc_i,k_crc3_i)
k_crc_i <- rbind(k_crc_i,k_crc4_i)

c198<-k_crc_i[,70] 
c364<-k_crc_i[,253] 
h198<-k_ctl1_i[,70] 
h364<-k_ctl1_i[,253] 
plot(c198,c364,xlim=c(0,1.3),ylim=c(0,2),col="red")
par(new = TRUE)
plot(h198, h364,xlim=c(0,1.3),ylim=c(0,2),col="green")

# boundary for capacity regions for Stage I, II, III, IV.
plot(k_crc1_i[,70],k_crc1_i[,253],xlim=c(0,3),ylim=c(0,3),col="blue",pch=1)
par(new = TRUE)
plot(k_crc2_i[,70],k_crc2_i[,253],xlim=c(0,3),ylim=c(0,3),col="orange",pch=2)
par(new = TRUE)
plot(k_crc3_i[,70],k_crc3_i[,253],xlim=c(0,3),ylim=c(0,3),col="red",pch=3)
par(new = TRUE)
#plot(k_crc4_i[,70],k_crc4_i[,253],xlim=c(0,7),ylim=c(0,7),col="brown",pch=4)
#par(new = TRUE)
plot(k_ctl1_i[,70],k_ctl1_i[,253],xlim=c(0,3),ylim=c(0,3),col="green",pch=5)
legend(2,3,legend=c("KLD for CRC SI samples", "KLD for CRC SII samples", "KLD for CRC SIII samples", "KLD for Healthy samples" ),col=c("blue","orange","red","green"),lty=1:1, cex=1.0)

plot(abs(colSums(k_crc_i,na.rm = T)-colSums(k_ctl1_i,na.rm = T)))
which(abs(colSums(k_crc_i,na.rm = T)-colSums(k_ctl1_i,na.rm = T))>22)

dim(k_ctl1_i)#43 700
hv<-rbind(ctl1[1:8,], ctl1[10:43,])
dim(k_crc_i)#  27 700
cv<-rbind(crc[1:21,], crc[23:27,])
names(hv) <- names(cv) 
hcTop20 <- rbind(cv, hv)
df<-scale(hcTop20)
df[is.nan(df)] <- 0
col <- colorRampPalette(brewer.pal(11, "RdYlBu"))(256)
hm <- heatmap(df, scale = "none", col =  col) 

dim(umi) # 26775   499
dim(umic)# 39865   499
dim(umii) #33320   499
umiB10<- array(0, dim = c(45*595*10,50))
umiiB10<- array(0, dim = c(56*595*10,50))
umicB10<-array(0,dim=c(68*595*10,50))
colB10<- array(0, dim = c(79*574*10,50))
ctl1B10<-array(0,dim=c(74*574*10,50))
for (i in 0:49){
  #i=2
  k=10*i+1
  if (i != 49) {
  umiB10[,i] <- umi[,k:(k+9)]
  umicB10[,i] <- umic[,k:(k+9)]
  umiiB10[,i] <- umii[,k:(k+9)]
  colB10[,i] <- colD22[,k:(k+9)]
  ctl1B10[,i] <- ctl1D22[,k:(k+9)]
  } else {
    umiB10[,i] <- umi[,(k-1):(k+8)]# repeat one FRl from the previous segment
    umicB10[,i] <- umic[,(k-1):(k+8)]# repeat one FRl from the previous segment
    umiiB10[,i] <- umii[,(k-1):(k+8)]# repeat one FRl from the previous segment
    colB10[,i] <- colD22[,(k-1):(k+8)]
    ctl1B10[,i] <- ctl1D22[,(k-1):(k+8)]
}
}





hgA <- hist(colB10[,19],breaks = 150 , plot = FALSE) # Save first histogram data
hgB <- hist(ctl1B10[,19] ,breaks = 150, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,500), ylim = c(0,14000)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,500), ylim = c(0,14000)) # Add 2nd histogram using different color





write.csv(colB10,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_col49_binned10frl.csv')
write.csv(ctl1B10,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_ctl1_74_binned10frl.csv')
write.csv(umiB10,'~/genomedk/matovanalysis/DELFI_analysis/python/umiseq_pon45_binned10frl.csv')
write.csv(umicB10,'~/genomedk/matovanalysis/DELFI_analysis/python/umiseq_CRpre68_binned10frl.csv')
write.csv(umiiB10,'~/genomedk/matovanalysis/DELFI_analysis/python/umiseq_pre56_binned10frl.csv')

umii356 <- umii[,356]
umii356 <-umii356[umii356>9]
 
umi356 <- umi[,356]
umi356 <-umi356[umi356>9]

write.csv(umi356,'~/genomedk/matovanalysis/DELFI_analysis/python/umiseq_pon45_frl356no0.csv')
write.csv(umii356,'~/genomedk/matovanalysis/DELFI_analysis/python/umiseq_pre56_frl356no0.csv')

hgA <- hist(umii356, breaks = 250 , plot = FALSE) # Save first histogram data
hgB <- hist(umi356,breaks = 150, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,180), ylim = c(0,1000)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,180), ylim = c(0,1000)) # Add 2nd histogram using different color



k_d2_col79_b10 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2_COL79_CTL1_binned10frl.csv')
kd2_col79_b10  <- k_d2_col79_b10[2:51,2]
plot(kd2_col79_b10)


k_umi_cPre68_b10 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceUMIseqCRPre68Pon45_binned10frl.csv')
k_umicPre68_b10  <- k_umi_cPre68_b10[2:51,2]
plot(k_umicPre68_b10)

k_umi_iPre56_b10 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceUMIseqPre56Pon45_binned10frl.csv')
k_umiiPre56_b10  <- k_umi_iPre56_b10[2:51,2]
plot(k_umiiPre56_b10)


# save all 57 PreOp I, convert 25 to bins (2 are very small)
# save all 68 PreOp C, convert 68 to bins
# save all 45 PON , convert 15 to bins
#samplesCTL3 <- sapply(m2$DELFI.ID[ctl3_list] , function(x) grep(x, x = pileupsC[9:90] )) 
listCRUK <- unlist(sapply(CRUKlist, function(x) grep(x, x = pileupsC [9:90])))
bams <-pileupsC[listCRUK ]
#bams <- list.files(pileupsC)
pileupsC_Bam <- list.files("~/genomedk/PolyA/faststorage/BACKUP/CRUK/plasma/N289", recursive = T, full.names = T, pattern = "_consensus.sort.bam$")
listCRUK_umi <- unlist(sapply(names(listCRUK), function(x) grep(x, x = pileupsC_Bam)))
bams <-pileupsC_Bam[listCRUK_umi ]
sapply(bams, function(x)file.copy(from=x, to = "~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/umiseq_paper/divergence/data/CRUK5Mb/"))

pileupsI_Bam <- list.files("~/genomedk/PolyA/faststorage/BACKUP/IMPROVE/sporacrc/pon/data/pon", recursive = T, full.names = T, pattern = "_consensus.sort.bam$")
sapply(pileupsI_Bam[1:46], function(x)file.copy(from=x, to = "~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/umiseq_paper/divergence/data/45PON5Mb/"))

It <- read.table("~/genomedk/matovanalysis/umiseq_analysis/IMPROVEptList",header = TRUE)
pileupsIM_Bam <- list.files("~/genomedk/PolyA/faststorage/BACKUP/IMPROVE/sporacrc/N227", recursive = T, full.names = T, pattern = "_consensus.sort.bam$")
It$index <- sapply(as.character(It$library_id), function(x) grep(x, pileupsIM_Bam)) 
preop_index <- unlist(It[ It$op_time_cat == -1, "index"]) # 56
pileupsIM_Bam[preop_index]
sapply(pileupsIM_Bam[preop_index], function(x)file.copy(from=x, to = "~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/umiseq_paper/divergence/data/57PRE5Mb/"))


       
       
       

pileupsD1 <- list.files("~/genomedk/DELFI1/Workspaces/maretty/frag_lens/5mb", recursive = T, full.names = T, pattern = "tsv")

listD1CRC27 <- sapply(iC, function(x) grep(x, x = pileupsD1 ))
pileupsD1[listD1CRC27] 
listD1CRC4_8 <- sapply(iC4, function(x) grep(x, x = pileupsD1 ))
pileupsD1[listD1CRC4_8] 
listD1CTL215 <- unlist(sapply(iH, function(x) grep(x, x = pileupsD1 )))




nbCRC <- length(iC)
crc <-matrix(,nrow=574,ncol=500)#nrow=555,ncol=702)
c <-matrix(,nrow=nbCRC ,ncol=499)#700)
nFc  <- vector()
# for all patients
for(i in 1:nbCRC) {
  #i= 1
  print(i)
    # replace with ~/DELFI1/Workspaces/maretty/frag_lens/ text filest
  auxFR <- read.table(pileupsD1[listD1CRC27[i]], header = TRUE) # sample per sample, file per file. 
  crc1 <- auxFR[,2:500]
    #crc1 <- bFr[as.character(bFr[,1])==as.character(iC[i]),]
    print(sum(crc1)) # 555 x 702
    
    if (i != 3){
    #nFc[i] <- sum(crc1[,3:702])
    }
    
    if (i == 10){
     # crcMax <- crc1[,3:702]
    }
    
    if (i == 16){
     # crcMin <- crc1[,3:702]
    }
    
    #c1 <- colSums(crc1[,3:702])
    #c [i,] <- c1
    # distribution for a bin crc1[1,3:702]
    if (i==1){
      crc <- crc1 
     # crcB <- t(crc1[,3:702])
    } else {
      names(crc) <- names(crc1) 
      crc <- rbind(crc , crc1)# agggregate merged results over all CRC samples
      if (i != 3){
     #   crcB <- rbind(crcB , t(crc1[,3:702]))# agggregate merged results over all CRC samples
        #ctl2 <- cbind(ctl , ctl1)# agggregate merged results over all CRC samples
      }
    }
}

k_crc1_i<-matrix(,nrow=7,ncol=700)
k_crc1_i3 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD1crc1_individual3.csv')
k_crc1_i[3,] <- k_crc1_i3[2:701,2]
plot(k_crc1_i[3,])
nbC1 <- length(iC1)
cs1 <-matrix(,nrow=555,ncol=702)
mC1=NULL
maC1=NULL
for(i in 1:nbC1) {
  #i= 1
  #cs11 <- bFr[bFr[,1]==iC1[i],]
  cs11 <- bFr[as.character(bFr[,1])==as.character(iC1[i]),]
  #sa_name <- paste0('~/genomedk/matovanalysis/DELFI_analysis/python/delfi1_crc1_individual', i,'.csv')
  #write.csv(cs11[,3:702],sa_name)
  #sa_name <- paste0('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD1crc1_individual', i,'.csv')
  #aux<- read.csv(sa_name)
  #k_crc1_i[i,] <- aux[2:701,2]
  #t1<-floor(i/9)
  #k <- (i-t1*9)
  #plot(k_crc1_i[i,],  ylim=range(c(0,5.2)), col= palette()[k], type="l") # there is something wrong with the colors
  #print(palette()[k])#print(col2rgb(2*i+50))
  #mC1[i]<-which.max(k_crc1_i[i,] )
  #maC1[i]<-max(k_crc1_i[i,] )
  #par(new = TRUE)
  # distribution for a bin crc1[1,3:702]
  if (i==1){
    cs1 <- cs11 
  } else {
    cs1 <- rbind(cs1 , cs11)# agggregate merged results over all CRC samples
  }
}
k_crc2_i<-matrix(,nrow=7,ncol=700)
k_crc2_i1 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD1crc2_individual1.csv')
k_crc2_i[1,] <- k_crc2_i1[2:701,2]
plot(k_crc2_i[1,])
nbC2 <- length(iC2)
cs2 <-matrix(,nrow=555,ncol=702)
mC1=NULL
maC1=NULL
# for all patients
for(i in 1:nbC2) {
  #i= 1
  #cs21 <- bFr[bFr[,1]==iC2[i],]
  cs21 <- bFr[as.character(bFr[,1])==as.character(iC2[i]),]
  sa_name <- paste0('~/genomedk/matovanalysis/DELFI_analysis/python/delfi1_crc2_individual', i,'.csv')
  #write.csv(cs21[,3:702],sa_name)
  sa_name <- paste0('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD1crc2_individual', i,'.csv')
  aux<- read.csv(sa_name)
  k_crc2_i[i,] <- aux[2:701,2]
  t1<-floor(i/9)
  k <- (i-t1*9)
  plot(k_crc2_i[i,],  ylim=range(c(0,5.5)), col= palette()[k], type="l") # there is something wrong with the colors
  print(palette()[k])#print(col2rgb(2*i+50))
  mC1[i]<-which.max(k_crc2_i[i,] )
  maC1[i]<-max(k_crc2_i[i,] )
  par(new = TRUE)
  # distribution for a bin crc1[1,3:702]
  if (i==1){
    cs2 <- cs21 
  } else {
    cs2 <- rbind(cs2 , cs21)# agggregate merged results over all CRC samples
  }
}


k_crc3_i<-matrix(,nrow=5,ncol=700)
k_crc3_i1 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD1crc3_individual1.csv')
k_crc3_i[1,] <- k_crc3_i1[2:701,2]
plot(k_crc3_i[1,])
mC1=NULL
maC1=NULL
nbC3 <- length(iC3)
cs3 <-matrix(,nrow=555,ncol=702)
# for all patients
for(i in 1:nbC3) {
  #i= 1
  #cs31 <- bFr[bFr[,1]==iC3[i],]
  cs31 <- bFr[as.character(bFr[,1])==as.character(iC3[i]),]
  sa_name <- paste0('~/genomedk/matovanalysis/DELFI_analysis/python/delfi1_crc3_individual', i,'.csv')
  #write.csv(cs31[,3:702],sa_name)
  sa_name <- paste0('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD1crc3_individual', i,'.csv')
  aux<- read.csv(sa_name)
  k_crc3_i[i,] <- aux[2:701,2]
  t1<-floor(i/9)
  k <- (i-t1*9)
  plot(k_crc3_i[i,],  ylim=range(c(0,6.2)), col= palette()[k], type="l") # there is something wrong with the colors
  print(palette()[k])#print(col2rgb(2*i+50))
  mC1[i]<-which.max(k_crc3_i[i,] )
  maC1[i]<-max(k_crc3_i[i,] )
  par(new = TRUE)
  # distribution for a bin crc1[1,3:702]
  if (i==1){
    cs3 <- cs31 
  } else {
    cs3 <- rbind(cs3 , cs31)# agggregate merged results over all CRC samples
  }
}
k_crc4_i<-matrix(,nrow=8,ncol=700)
k_crc4_i1 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD1crc4_individual1.csv')
k_crc4_i[1,] <- k_crc4_i1[2:701,2]
plot(k_crc4_i[1,])
mC1=NULL
maC1=NULL
nbC4 <- length(iC4)
cs4 <-matrix(,nrow=555,ncol=702)
cs4 <-matrix(,nrow=574,ncol=499)
# for all patients
for(i in 1:nbC4) {
  #i= 1
  #cs41 <- bFr[bFr[,1]==iC4[i],]
  #cs41 <- bFr[as.character(bFr[,1])==as.character(iC4[i]),]
  auxFR <- read.table(pileupsD1[listD1CRC4_8[i]], header = TRUE) # sample per sample, file per file. 
  cs41 <- auxFR[,2:500]
  sa_name <- paste0('~/genomedk/matovanalysis/DELFI_analysis/python/delfi1_crc4_individual', i,'.csv')
  #write.csv(cs41[,3:702],sa_name)
  if (i != 3){
    #sa_name <- paste0('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD1crc4_individual', i,'.csv')
    #aux<- read.csv(sa_name)
    #k_crc4_i[i,] <- aux[2:701,2]
    #t1<-floor(i/9)
    #k <- (i-t1*9)
    #plot(k_crc4_i[i,],  ylim=range(c(0,7.2)), col= palette()[k], type="l") # there is something wrong with the colors
    #print(palette()[k])#print(col2rgb(2*i+50))
    #mC1[i]<-which.max(k_crc4_i[i,] )
    #maC1[i]<-max(k_crc4_i[i,] )
    #par(new = TRUE)
  }
  # distribution for a bin crc1[1,3:702]
  if (i==1){
    cs4 <- cs41 
  } else {
    names(cs4)<-names(cs41)
    cs4 <- rbind(cs4 , cs41)# agggregate merged results over all CRC samples
  }
}

cS <-c("PGDX5881P"  ,"PGDX5882P",  "PGDX5883P1", "PGDX5884P"  ,"PGDX5886P"  ,"PGDX5889P"  ,"PGDX5890P1" ,"PGDX5891P"  ,"PGDX5892P",  "PGDX5894P1",
       "PGDX5969P1" ,"PGDX5895P1",
"PGDX5899P1", "PGDX5900P1" ,"PGDX5903P1", "PGDX5963P" , "PGDX5964P" , "PGDX5965P" , "PGDX5966P",  "PGDX5967P"  ,"PGDX5968P",
"PGDX6249P1", "PGDX8820P1",
"PGDX8823P" , "PGDX8825P1", "PGDX8828P",  "PGDX8829P1")
dimnames(c)[[1]]<-cS
cF <- c(191, 192 ,194 ,195, 196, 197, 198, 200, 201, 203 ,204, 205 ,206 ,357, 358 ,362 ,364 ,367, 368 ,370)
c20<- c[,cL]
dimnames(c20)[[2]]<-cF

k_ctl1_i<-matrix(,nrow=43,ncol=700)
k_ctl1_i3 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD1ctl1_individual3.csv')
k_ctl1_i[3,] <- k_ctl1_i3[2:701,2]
plot(k_ctl1_i[3,])

nbCTL <- length(iH)#/5#*0.4#/ 5 # devide by 5 until laptop comes
#ctl <-matrix(,nrow=555,ncol=702)
ctlB <-matrix(,nrow=700,ncol=555)
h <-matrix(,nrow=nbCTL ,ncol=700)
nFh  <- vector()
#nFh <- c(nFh, 1:nbCTL)
# for all patients
mC1=NULL
maC1=NULL
totFRs_D1_CTL1<- vector()
totFRs_D1_CTL11<- vector()


ctl <-matrix(,nrow=574,ncol=500)#nrow=555,ncol=702)
c <-matrix(,nrow=nbCTL ,ncol=499)#700)

for(i in 1:nbCTL) {#nbCTL
 # i= 25
  print(i)
  auxFR <- read.table(pileupsD1[unlist(listD1CTL215)[i]], header = TRUE) # sample per sample, file per file. 
  ctl1 <- auxFR[,2:500]
 #ctl1 <- bFr[as.character(bFr[,1])==as.character(iH[i]),]
  print(sum(ctl1))
  #sa_name <- paste0('~/genomedk/matovanalysis/DELFI_analysis/python/delfi1_ctl1_individual', i,'.csv')
  #write.csv(ctl1[,3:702],sa_name)
  if (i != 9){
    #totFRs_D1_CTL1[i] <- sum(ctl1[1:499,3:702])
    #print(totFRs_D1_CTL1[i])
    #totFRs_D1_CTL11[i] <- sum(ctl1[,3:702])
    
    #sa_name <- paste0('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD1ctl1_individual', i,'.csv')
    #aux<- read.csv(sa_name)
    #k_ctl1_i[i,] <- aux[2:701,2]
    #t1<-floor(i/9)
    #k <- (i-t1*9)
    #plot(k_ctl1_i[i,],  ylim=range(c(0,2.8)), col= palette()[k], type="l") # there is something wrong with the colors
    #print(palette()[k])#print(col2rgb(2*i+50))
    #mC1[i]<-which.max(k_ctl1_i[i,] )
    #maC1[i]<-max(k_ctl1_i[i,] )
   #par(new = TRUE)
  }
  if (i == 10){
    #ctlMax <- ctl1[,3:702]
  }
  
  if (i == 11){
    #ctlMin <- ctl1[,3:702]
  }
  
  
  if (i != 9){
   # nFh[i] <- sum(ctl1[,3:702])
  }
  
  #h1 <- colSums(ctl1[,3:702])
  #h[i,] <- h1
  # distribution for a bin ctl1[1,3:702]
  if (i==1){
    ctl <- ctl1 
   # ctlB <- t(ctl1[,3:702])
  } else {
    names(ctl)<-names(ctl1)
    ctl <- rbind(ctl , ctl1)# agggregate merged results over all CRC samples
    if (i != 9){
    #ctlB <- rbind(ctlB , t(ctl1[,3:702]))# agggregate merged results over all CRC samples
    #ctl2 <- cbind(ctl , ctl1)# agggregate merged results over all CRC samples
    }
  }
}
write.csv(ctl,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi1_ctl86HG38.csv')
write.csv(ctl,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi1_ctl129HG38.csv')
write.csv(crc,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi1_crc27HG38.csv')
write.csv(ctl,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi1_ctl43HG38.csv')
k_cHG38 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceCRC_HG38.csv')
kcHG38<-k_cHG38[2:500,2]
plot(kcHG38)
which.max(kcHG38)

for(i in (nbCTL+1):(nbCTL+nbCTL+1)) {
  #i= 1
  ctl1 <- bFr[bFr[,1]==iH[i],]
  
  print(dim(ctl1))
  
  # distribution for a bin ctl1[1,3:702]
  if (i==(nbCTL+1)){
    ctl <- ctl1 
  } else {
    ctl <- rbind(ctl , ctl1)# agggregate merged results over all CRC samples
  }
}
#compare all CTL vs all CRC samples in terms of overall number of fragments, i.e. coverage
nFh[9]<-nFh[43]
nFh<-nFh[1:42]
nFc[3]<-nFc[27]
nFc<-nFc[1:26]
plot(nFh, col="green",ylim=range(c(23000000,86000000)))
par(new = TRUE)
plot(nFc, col="red",ylim=range(c(23000000,86000000)))

#hS <-c("PGDX5881P"  ,"PGDX5882P",  "PGDX5883P1", "PGDX5884P"  ,"PGDX5886P"  ,"PGDX5889P"  ,"PGDX5890P1" ,"PGDX5891P"  ,"PGDX5892P",  "PGDX5894P1",
       #"PGDX5969P1" ,"PGDX5895P1",
       #"PGDX5899P1", "PGDX5900P1" ,"PGDX5903P1", "PGDX5963P" , "PGDX5964P" , "PGDX5965P" , "PGDX5966P",  "PGDX5967P"  ,"PGDX5968P",
       #"PGDX6249P1", "PGDX8820P1",
       #"PGDX8823P" , "PGDX8825P1", "PGDX8828P",  "PGDX8829P1")
dimnames(h)[[1]]<-hS
hF <- c(191, 192 ,194 ,195, 196, 197, 198, 200, 201, 203 ,204, 205 ,206 ,357, 358 ,362 ,364 ,367, 368 ,370)
h20<- h[,cL]
dimnames(h20)[[2]]<-hF

# for all bins
# extract FR histograms and assign to a small matrix
# repeat with Healthy
# plug two matrices into KLD


length(unique(sName))# 533 patients, 295815 bins in total, which is 555 bin per patient, which covers 2.775G of genome

ln <- as.integer(unlist(lcc[,3:702]))# convert to numbers
lnm <- matrix(ln, ncol = dim(lcc)[1], byrow = (dim(lcc)[2]-2)) # convert back to matrix form

dn <- as.integer(unlist(dcc[,3:702]))# convert to numbers
dnm <- matrix(dn, ncol = dim(dcc)[1], byrow = (dim(dcc)[2]-2)) # convert back to matrix form

bn <- as.integer(unlist(bcc[,3:702]))# convert to numbers
bnm <- matrix(bn, ncol = dim(bcc)[1], byrow = (dim(bcc)[2]-2)) # convert back to matrix form

gn <- as.integer(unlist(gcc[,3:702]))# convert to numbers
gnm <- matrix(gn, ncol = dim(gcc)[1], byrow = (dim(gcc)[2]-2)) # convert back to matrix form

pn <- as.integer(unlist(pcc[,3:702]))# convert to numbers
pnm <- matrix(pn, ncol = dim(pcc)[1], byrow = (dim(pcc)[2]-2)) # convert back to matrix form

on <- as.integer(unlist(ovc[,3:702]))# convert to numbers
onm <- matrix(on, ncol = dim(ovc)[1], byrow = (dim(ovc)[2]-2)) # convert back to matrix form

cn <- as.integer(unlist(crc))#[,3:702]))# convert to numbers
cnm <- matrix(cn, ncol = dim(crc)[1], byrow = (dim(crc)[2]-2)) # convert back to matrix form
cnmD<-cnm*2 # mimic double coverage
cnmH<-cnm/2 # mimic half coverage

cn1 <- as.integer(unlist(cs1[,3:702]))# convert to numbers
cnm1 <- matrix(cn1, ncol = dim(cs1)[1], byrow = (dim(cs1)[2]-2)) # convert back to matrix form

cn2 <- as.integer(unlist(cs2[,3:702]))# convert to numbers
cnm2 <- matrix(cn2, ncol = dim(cs2)[1], byrow = (dim(cs2)[2]-2)) # convert back to matrix form

cn3 <- as.integer(unlist(cs3[,3:702]))# convert to numbers
cnm3 <- matrix(cn3, ncol = dim(cs3)[1], byrow = (dim(cs3)[2]-2)) # convert back to matrix form

cn4 <- as.integer(unlist(cs4[,3:702]))# convert to numbers
cnm4 <- matrix(cn4, ncol = dim(cs4)[1], byrow = (dim(cs4)[2]-2)) # convert back to matrix form

hn <- as.integer(unlist(ctl[,3:702]))# convert to numbers
hnm <- matrix(hn, ncol = dim(ctl)[1], byrow = (dim(ctl)[2]-2)) # convert back to matrix form

#compare all bins CTL vs all bins CRC per FRL
plot(rowSums(hnm), col="green",ylim=range(c(0,55000000)))
par(new = TRUE)
plot(rowSums(cnm), col="red",ylim=range(c(0,55000000)))

#KL.divergence(t(cnm), t(hnm), k = 10, algorithm=c("kd_tree", "cover_tree", "brute")) # from FNN; input must be matrices
# 504.4316 507.1018 508.6900 509.2213 509.0071 508.4851 507.8492 507.7506 507.1370 506.2944

#KLD(t(cnm), t(hnm)) # from LaplacesDemon
write.csv(t(v364[sample(1:14430,14430)]),'~/genomedk/matovanalysis/DELFI_analysis/python/delfi1_c364_14430.csv')
write.csv(t(hnm[364,]),'~/genomedk/matovanalysis/DELFI_analysis/python/delfi1_h364_23310.csv')

write.csv(t(lnm),'G:\\matovanalysis/DELFI_analysis/python/delfi1_lcc.csv')
write.csv(t(dnm),'G:\\matovanalysis/DELFI_analysis/python/delfi1_dcc.csv')
write.csv(t(bnm),'G:\\matovanalysis/DELFI_analysis/python/delfi1_bcc.csv')
write.csv(t(gnm),'G:\\matovanalysis/DELFI_analysis/python/delfi1_gcc.csv')
write.csv(t(pnm),'G:\\matovanalysis/DELFI_analysis/python/delfi1_pcc.csv')
write.csv(t(onm),'G:\\matovanalysis/DELFI_analysis/python/delfi1_ovc.csv')
write.csv(t(cnm),'G:\\matovanalysis/DELFI_analysis/python/delfi1_crc.csv')
write.csv(t(cnmD),'G:\\matovanalysis/DELFI_analysis/python/delfi1_crcD.csv')# double
write.csv(t(cnmH),'G:\\matovanalysis/DELFI_analysis/python/delfi1_crcH.csv')# half- doesnt work well
write.csv(crcMax,'G:\\matovanalysis/DELFI_analysis/python/delfi1_crcMax.csv')# max coverage crc pt
write.csv(crcMin,'G:\\matovanalysis/DELFI_analysis/python/delfi1_crcMin.csv')# min coverage crc pt
write.csv(ctlMax,'G:\\matovanalysis/DELFI_analysis/python/delfi1_ctlMax.csv')# max coverage crc pt
write.csv(ctlMin,'G:\\matovanalysis/DELFI_analysis/python/delfi1_ctlMin.csv')# min coverage crc pt
write.csv(t(hnm),'G:\\matovanalysis/DELFI_analysis/python/delfi1_ctl.csv')
write.csv(t(hnm),'G:\\matovanalysis/DELFI_analysis/python/delfi1_ctl2.csv')
ctlBB<-unname(ctlB)
write.csv(ctlBB,'G:\\matovanalysis/DELFI_analysis/python/delfi1_ctlB.csv')
crcBB<-unname(crcB)
write.csv(crcBB,'G:\\matovanalysis/DELFI_analysis/python/delfi1_crcB.csv')
write.csv(hM364,'G:\\matovanalysis/DELFI_analysis/python/delfi1_ctl364.csv')
write.csv(cM364,'G:\\matovanalysis/DELFI_analysis/python/delfi1_crc364.csv')
write.csv(hM205,'G:\\matovanalysis/DELFI_analysis/python/delfi1_ctl205.csv')
write.csv(cM205,'G:\\matovanalysis/DELFI_analysis/python/delfi1_crc205.csv')
write.csv(hM198,'G:\\matovanalysis/DELFI_analysis/python/delfi1_ctl198.csv')
write.csv(cM198,'G:\\matovanalysis/DELFI_analysis/python/delfi1_crc198.csv')
write.csv(t(cnm1),'G:\\matovanalysis/DELFI_analysis/python/delfi1_crc1.csv')
write.csv(t(cnm2),'G:\\matovanalysis/DELFI_analysis/python/delfi1_crc2.csv')
write.csv(t(cnm3),'G:\\matovanalysis/DELFI_analysis/python/delfi1_crc3.csv')
write.csv(t(cnm4),'G:\\matovanalysis/DELFI_analysis/python/delfi1_crc4.csv')

write.csv(t(v364[sample(1:14430,100)]),'~/genomedk/matovanalysis/DELFI_analysis/python/delfi1_v364_100.csv')

write.csv(he1_bins,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi1_ctl1Bins.csv')
write.csv(crc_bins,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi1_crcBins.csv')
k_cBins <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceCRC_Bins.csv')
kcBins<-k_cBins[2:701,2]
plot(kcBins)



k_ch2 <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRCh2.csv')
k_cmaxmin <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRCmaxmin.csv')

kcmax_hmax <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCTLmaxCRCmax.csv')

k_hmax <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCTLmax.csv')
k_hmin <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCTLmin.csv')

k_cbin2 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceCRC_bins364.csv')

k_cbin <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceCRCbins.csv')
k_cmax <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRCmax.csv')
k_cmin <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRCmin.csv')
k_c12 <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC1CRC2.csv')
k_c13 <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC1CRC3.csv')
k_c14 <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC1CRC4.csv')
k_c23 <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC2CRC3.csv')
k_c24 <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC2CRC4.csv')
k_c34 <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC3CRC4.csv')
k_cD <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRCd.csv')
k_cH <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRCh.csv')
k_c1 <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC1.csv')
k_c2 <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC2.csv')
k_c3 <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC3.csv')
k_c4 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceCRC4.csv')
k_c364 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceCRC364.csv')
k_c205 <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC205.csv')
k_c198 <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC198.csv')
k_l <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceLCC.csv')
k_cl <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC_LCC.csv')
k_d <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceDCC.csv')
k_cd <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC_DCC.csv')
k_b <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceBCC.csv')
k_cb <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC_BCC.csv')
k_cg <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC_GCC.csv')
k_g <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceGCC.csv')
k_p <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergencePCC.csv')
k_c <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceCRC.csv') # CRC
kc<-k_c[2:701,2]
plot(kc)
k_cR <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceCRC_Reverse.csv') # CRC  REVERSE

k_o <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceOVC.csv')
k_co <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC_OVC.csv')
k_cp <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC_PCC.csv')

kcbin2<-k_cbin2[2:701,2]
plot(kcbin2)
kcmaxmin<-k_cmaxmin[2:701,2]
plot(kcmaxmin)
kch2<-k_ch2[2:701,2]
plot(kch2)
khmaxcmax<-kcmax_hmax[2:701,2]
plot(khmaxcmax)
khmin<-k_hmin[2:701,2]
plot(khmin)
khmax<-k_hmax[2:701,2]
plot(khmax)
#kcbin<-k_cbin[2:701,2]
plot(k_cbin)
kcmax<-k_cmax[2:701,2]
plot(kcmax)
kcmin<-k_cmin[2:701,2]
plot(kcmin)
kcD<-k_cD[2:701,2]
plot(kcD)
kc23<-k_c23[2:701,2]
plot(kc23)
kc24<-k_c24[2:701,2]
plot(kc24)
kc34<-k_c34[2:701,2]
plot(kc34)
kc14<-k_c14[2:701,2]
plot(kc14)
kc13<-k_c13[2:701,2]
plot(kc13)
kc12<-k_c12[2:701,2]
plot(kc12)
kc1<-k_c1[2:701,2]
plot(kc1)
kc2<-k_c2[2:701,2]
plot(kc2)
kc3<-k_c3[2:701,2]
plot(kc3)
kc4<-k_c4[2:701,2]
plot(kc4)
kc364<-k_c364[2:701,2]
plot(kc364)
kc205<-k_c205[2:701,2]
plot(kc205)
kc198<-k_c198[2:701,2]
plot(kc198)
kl<-k_l[2:701,2]
plot(kl)
kcl<-k_cl[2:701,2]
plot(kcl)
kd<-k_d[2:701,2]
plot(kd)
kcd<-k_cd[2:701,2]
plot(kcd)
kb<-k_b[2:701,2]
plot(kb)
kcb<-k_cb[2:701,2]
plot(kcb)
kcg<-k_cg[2:701,2]
plot(kcg)
kg<-k_g[2:701,2]
plot(kg)
kcp<-k_cp[2:701,2]
plot(kcp)
kp<-k_p[2:701,2]
plot(kp)
kc<-k_c[2:701,2]
plot(kc)
kcr<-k_cR[2:701,2]
plot(kcr)
colorKLD <- rowSums(hnm)/rowSums(cnm)
div <- data.frame(kc) 
div$clr = colorKLD
ve <- seq(1, 700, by=1)
cK <- ggplot(div, aes(x= ve, y = kc, color = clr) )  + geom_point()
cK+scale_color_gradientn(colours = rainbow(5))
ko<-k_o[2:701,2]
plot(ko)
kco<-k_co[2:701,2]
plot(kco)

plot(khmin, ylim=range(c(0,4)), col="green", main = "Divergence [in bits] from Healthy of healthy with highest and lowest coverage", type="l")
par(new = TRUE)
plot(khmax, ylim=range(c(0,4)), col="red", type="l")
legend(500,3,legend=c("CTL pt w 64mln fragments", "CTL w 28mln fragments"),col=c("green","red"),lty=1:1, cex=1.0)

plot(kcmin, ylim=range(c(0,4)), col="orange", main = "Divergence [in bits] from Healthy of individual samples", type="l")
par(new = TRUE)
plot(kcmax, ylim=range(c(0,4)), col="red", type="l")
par(new = TRUE)
plot(khmin, ylim=range(c(0,4)), col="green", type="l")
par(new = TRUE)
plot(khmax, ylim=range(c(0,4)), col="blue", type="l")
legend(420,3,legend=c("CRC pt w 23mln fragments", "CRC pt w 85mln fragments","CTL pt w 28mln fragments", "CTL w 64mln fragments"),col=c("orange","red","green","blue"),lty=1:1, cex=1.0)

plot(kcmin, ylim=range(c(0,4)), col="green", main = "Divergence [in bits] from Healthy of individual samples", type="l")
par(new = TRUE)
plot(kcmax, ylim=range(c(0,4)), col="red", type="l")
legend(500,3,legend=c("CRC pt w 23mln fragments", "CRC pt w 85mln fragments"),col=c("green","red"),lty=1:1, cex=1.0)

plot(kc23, ylim=range(c(0,2.4)), col="green", main = "Divergence [in bits] from CRC II of CRC III,IV per fragment length [in bp]", type="l")
par(new = TRUE)
plot(kc24, ylim=range(c(0,2.4)), col="blue", type="l")
par(new = TRUE)
plot(kc34, ylim=range(c(0,2.4)), col="red", type="l")
legend(500,1.6,legend=c("CRC Stage II->III", "CRC Stage II->IV", "CRC Stage III->IV"),col=c("green","blue","red"),lty=1:1, cex=1.0)

plot(kc12, ylim=range(c(0,2.2)), col="green", main = "Divergence [in bits] from CRC I of CRC II-IV per fragment length [in bp]", type="l")
par(new = TRUE)
plot(kc13, ylim=range(c(0,2.2)), col="blue", type="l")
par(new = TRUE)
plot(kc14, ylim=range(c(0,2.2)), col="red", type="l")
legend(500,1.6,legend=c("CRC Stage I->II", "CRC Stage I->III", "CRC Stage I->IV"),col=c("green","blue","red"),lty=1:1, cex=1.0)

plot(kc1, ylim=range(c(0,5)), col="green", main = "Divergence [in bits] from healthy of CRC Stages I-IV per fragment length [in bp]", type="l")
par(new = TRUE)
plot(kc2, ylim=range(c(0,5)), col="blue", type="l")
par(new = TRUE)
plot(kc3, ylim=range(c(0,5)), col="orange", type="l")
par(new = TRUE)
plot(kc4, ylim=range(c(0,5)), col="red", type="l")
legend(500,2.6,legend=c("CRC Stage I", "CRC Stage II", "CRC Stage III", "CRC Stage IV"),col=c("green","blue","orange","red"),lty=1:1, cex=1.0)

plot(kc, ylim=range(c(0,2.7)), col="green", main = "Divergence [in bits] from healthy for cancers per fragment length [in bp]", type="l")
par(new = TRUE)
plot(ko, ylim=range(c(0,2.7)), col="blue", type="l")
par(new = TRUE)
plot(kb, ylim=range(c(0,2.7)), col="orange", type="p")
par(new = TRUE)
plot(kp, ylim=range(c(0,2.7)), col="red", type="l")
par(new = TRUE)
plot(kg, ylim=range(c(0,2.7)), col="brown", type="l")
par(new = TRUE)
plot(kl, ylim=range(c(0,2.7)), col="pink", type="l")
par(new = TRUE)
plot(kd, ylim=range(c(0,2.7)), col="purple", type="l")
legend(500,2.6,legend=c("CRC", "Ovarian", "Breast", "Pancreatic", "Gastric", "Lung", "Bile duct"), col=c("green","blue","orange","red","brown","pink","purple"),lty=1:1, cex=1.0)

plot(kcd, ylim=range(c(0,1.7)), col="green", main = "Divergence [in bits] from CRC of other cancers per fragment length [in bp]", type="l")
par(new = TRUE)
plot(kco, ylim=range(c(0,1.7)), col="blue", type="l")
par(new = TRUE)
plot(kcb, ylim=range(c(0,1.7)), col="orange", type="p")
par(new = TRUE)
plot(kcp, ylim=range(c(0,1.7)), col="red", type="l")
par(new = TRUE)
plot(kcg, ylim=range(c(0,1.7)), col="brown", type="l")
par(new = TRUE)
plot(kcl, ylim=range(c(0,1.7)), col="purple", type="l")
legend(500,1.7,legend=c("Bile duct", "Ovarian", "Breast", "Pancreatic", "Gastric", "Lung"), col=c("green","blue","orange","red","brown","purple"),lty=1:1, cex=1.0)

sum(s1[which(kc>1)])/sum(s1) # % fragments belonging to KLD>1 based on one bin analysis only

sum(cnm[which(kc>2.8),])/sum(cnm) # 0.0001387829
sum(cnm[which(kc>2),])/sum(cnm) # 0.06610582
sum(cnm[which(kc>1),])/sum(cnm) # 0.1513125

# compare Healthy and CRC number of fragments per bin for FR=364 and 222
hist(cnm[364,], xlim = c(0,70), breaks = 20)
hist(hnm[364,], xlim = c(0,70), breaks = 20)

# 1st most discriminate bin for 205 bp
hgA <- hist(cM205[,369],xlim = c(0,70), breaks = 20 , plot = FALSE) # Save first histogram data
hgB <- hist(hM205[,369], xlim = c(0,70),breaks = 20, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,80), ylim = c(0,4)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,80), ylim = c(0,4)) # Add 2nd histogram using different color

# 1st most discriminate bin for 198 bp
hgA <- hist(cM198[,447], breaks = 10 , plot = FALSE) # Save first histogram data
hgB <- hist(hM198[,447], breaks = 25, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,600), ylim = c(0,9)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,600), ylim = c(0,9)) # Add 2nd histogram using different color

# 1st most discriminate bin for 364 bp
hgA <- hist(cM364[,308],xlim = c(0,70), breaks = 20 , plot = FALSE) # Save first histogram data
hgB <- hist(hM364[,308], xlim = c(0,70),breaks = 50, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,40), ylim = c(0,6)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,40), ylim = c(0,6)) # Add 2nd histogram using different color

# 2nd most discriminate bin for 364 bp
hgA <- hist(cM364[,519],xlim = c(0,70), breaks = 10 , plot = FALSE) # Save first histogram data
hgB <- hist(hM364[,519], xlim = c(0,70),breaks = 50, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,50), ylim = c(0,6)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,50), ylim = c(0,6)) # Add 2nd histogram using different color

hgA <- hist(ctlMin[,364],xlim = c(0,70), breaks = 20 , plot = FALSE) # Save first histogram data
hgB <- hist(hnm[364,], xlim = c(0,70),breaks = 50, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,70), ylim = c(0,2200)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,70), ylim = c(0,2200)) # Add 2nd histogram using different color

hgA <- hist(crcMin[,364],xlim = c(0,70), breaks = 20 , plot = FALSE) # Save first histogram data
hgB <- hist(hnm[364,], xlim = c(0,70),breaks = 50, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,70), ylim = c(0,2200)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,70), ylim = c(0,2200)) # Add 2nd histogram using different color

hgA <- hist(cnm[364,],xlim = c(0,70), breaks = 50, plot = FALSE) # Save first histogram data
hgB <- hist(hnm[364,], xlim = c(0,70),breaks = 50, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,70), ylim = c(0,2200)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,70), ylim = c(0,2200)) # Add 2nd histogram using different color

hgA <- hist(cnm[137,],xlim = c(0,70), breaks = 100, plot = FALSE) # Save first histogram data
hgB <- hist(hnm[137,], xlim = c(0,70),breaks = 50, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,1200), ylim = c(0,5000)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,1200), ylim = c(0,5000)) # Add 2nd histogram using different color

hgA <- hist(crc[,137],xlim = c(0,70), breaks = 150, plot = FALSE) # Save first histogram data
hgB <- hist(ctl[,137], xlim = c(0,70),breaks = 50, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,1200), ylim = c(0,2500)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,1200), ylim = c(0,2500)) # Add 2nd histogram using different color

hgA <- hist(cs4[,137],xlim = c(0,70), breaks = 150, plot = FALSE) # Save first histogram data
hgB <- hist(ctl[,137], xlim = c(0,70),breaks = 50, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,1200), ylim = c(0,2500)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,1200), ylim = c(0,2500)) # Add 2nd histogram using different color

hgA <- hist(cnm4[137,],xlim = c(0,70), breaks = 100, plot = FALSE) # Save first histogram data
hgB <- hist(hnm[137,], xlim = c(0,70),breaks = 50, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,1200), ylim = c(0,5000)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,1200), ylim = c(0,5000)) # Add 2nd histogram using different color

hgA <- hist(cnm[277,],xlim = c(0,70), breaks = 100, plot = FALSE) # Save first histogram data
hgB <- hist(hnm[277,], xlim = c(0,70),breaks = 50, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,70), ylim = c(0,5000)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,70), ylim = c(0,5000)) # Add 2nd histogram using different color

hgA <- hist(crc[,277],xlim = c(0,70), breaks = 150, plot = FALSE) # Save first histogram data
hgB <- hist(ctl[,277], xlim = c(0,70),breaks = 40, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,75), ylim = c(0,2500)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,75), ylim = c(0,2500)) # Add 2nd histogram using different color

hgA <- hist(cs4[,277],xlim = c(0,70), breaks = 150, plot = FALSE) # Save first histogram data
hgB <- hist(ctl[,277], xlim = c(0,70),breaks = 40, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,75), ylim = c(0,2500)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,75), ylim = c(0,2500)) # Add 2nd histogram using different color

hgA <- hist(cnm4[277,],xlim = c(0,70), breaks = 100, plot = FALSE) # Save first histogram data
hgB <- hist(hnm[277,], xlim = c(0,70),breaks = 50, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,80), ylim = c(0,5000)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,80), ylim = c(0,5000)) # Add 2nd histogram using different color


hgA <- hist(cnm3[210,],xlim = c(0,70), breaks = 50, plot = FALSE) # Save first histogram data
hgB <- hist(hnm[210,], xlim = c(0,70),breaks = 50, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,70), ylim = c(0,2200)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,70), ylim = c(0,2200)) # Add 2nd histogram using different color

hgA <- hist(cnm[205,],breaks = 30, plot = FALSE) # Save first histogram data
hgB <- hist(hnm[205,], breaks = 60, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,560),  ylim = c(0,3800)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,560),  ylim = c(0,3800)) # Add 2nd histogram using different color

hgA <- hist(cnm[168,],xlim = c(0,15000), breaks = 90, plot = FALSE) # Save first histogram data
hgB <- hist(hnm[168,], xlim = c(0,15000),breaks = 60, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,4000), ylim = c(0,3200)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,4000), ylim = c(0,3200)) # Add 2nd histogram using different color

hgA <- hist(cnm[164,],xlim = c(0,15000), breaks = 90, plot = FALSE) # Save first histogram data
hgB <- hist(hnm[164,], xlim = c(0,15000),breaks = 60, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,4200), ylim = c(0,3200)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,4200), ylim = c(0,3200)) # Add 2nd histogram using different color

hist(hnm[168,])
hist(cnm[168,])

#SHOW A HEATMAP OF HOW NBS FOR 364 vary per bin (555bins) per CRC patient or normal
# 555x27 and 555x43
cM364 <- matrix(cnm[364,], ncol = 555, byrow = 26) #  
hM364 <- matrix(hnm[364,], ncol = 555, byrow = 42) # 
cM205 <- matrix(cnm[205,], ncol = 555, byrow = 26) #  
hM205 <- matrix(hnm[205,], ncol = 555, byrow = 42) # 
cM198 <- matrix(cnm[198,], ncol = 555, byrow = 26) #  
hM198 <- matrix(hnm[198,], ncol = 555, byrow = 42) # 

# compute divergence over all bins
crc_bins <- matrix(cnm, ncol = 555, byrow = 26) #  
he1_bins <- matrix(hnm, ncol = 555, byrow = 42) # 

# check variance across bins/genomic positions to find highest variability/difference CRC to CTL
vC364<-apply(cM364,2,var)
vH364<-apply(hM364,2,var)
mC364<-apply(cM364,2,mean)
mH364<-apply(hM364,2,mean)

plot(vC364, ylim=range(c(0,200)), col="red")
par(new = TRUE)
plot(vH364, ylim=range(c(0,200)), col="green")
plot(mC364, ylim=range(c(0,30)), col="red")
par(new = TRUE)
plot(mH364, ylim=range(c(0,30)), col="green")

dV364 <- vH364-vC364
dM364 <- mH364-mC364
plot(dV364)
which(dV364==max(dV364))# for FRL364, CRL is more dispersed in bin #475
which(dV364==min(dV364)) # for FRL364, CRC is more dispersed in bin #146
rV364 <- vH364/vC364
which(rV364>8) # bins # 74 220 282 373 444 453 455 495 498 525 553
which(rV364>9) # bins # 74 220 453 495
plot(dM364)
which(dM364==max(dM364))# for FRL364, CRL is more dispersed in bin #475
which(dM364==min(dM364)) # for FRL364, CRC is more dispersed in bin #146

#SHOW A HEATMAP OF HOW NBS FOR 202 vary per bin (555bins) per CRC patient or normal
# 555x27 and 555x43
cM205 <- matrix(cnm[205,], ncol = 555, byrow = 26) #  
hM205 <- matrix(hnm[205,], ncol = 555, byrow = 42) # 
# check variance across bins/genomic positions to find highest variability/difference CRC to CTL
vC205<-apply(cM205,2,var)
vH205<-apply(hM205,2,var)

plot(vC205, ylim=range(c(0,5000)), col="red")
par(new = TRUE)
plot(vH205, ylim=range(c(0,5000)), col="green")

dV205 <- vH205-vC205
plot(dV205)
which(dV205==max(dV205))# for FRL364, CRL is more dispersed in bin #475
which(dV205==min(dV205)) # for FRL364, CRC is more dispersed in bin #146
rV205 <- vH205/vC205

cM198 <- matrix(cnm[198,], ncol = 555, byrow = 26) #  
hM198 <- matrix(hnm[198,], ncol = 555, byrow = 42) # 
# check variance across bins/genomic positions to find highest variability/difference CRC to CTL
vC198<-apply(cM198,2,var)
vH198<-apply(hM198,2,var)

plot(vC198, ylim=range(c(0,12000)), col="red")
par(new = TRUE)
plot(vH198, ylim=range(c(0,12000)), col="green")

dV198 <- vH198-vC198
plot(dV198)
which(dV198==max(dV198))# for FRL364, CRL is more dispersed in bin #475
which(dV198==min(dV198)) # for FRL364, CRC is more dispersed in bin #146

# mix for clustering test
c364 <- cnm[364,]
h364 <- hnm[364,]

ze36 = c (1,  16, 14,  9,   12,  4,  4,  6,  7,   14, 15 , 9, 2, 2, 3, 5, 6, 7, 8, 9, 9, 11, 14, 19, 19 ,2, 4, 11, 14, 9, 14, 1, 1,
        1, 1 ,1 ,2 ,2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10,
        10, 10, 11, 11, 12, 12, 14, 14, 13, 14, 15, 16, 1, 16, 17, 17, 18, 18, 1, 2, 2, 21, 4, 9, 11, 1, 1, 1, 1, 1, 1, 1, 2, 2,2,  
         3, 3, 3, 3, 3, 3, 4, 4,  4, 4, 4, 4, 4, 4 ,4, 5, 5, 6, 6, 6, 6,6, 7,  7, 7, 7, 8, 8, 8,9, 9, 9, 9, 10, 10, 10, 10, 10, 10,
        10, 11, 12, 12, 12, 12, 12,13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 15, 15, 15, 17, 17, 17, 18, 19, 19, 20, 20, 21, 21, 9)
 
cL <- which(kc>2.38) # 191 192 194 195 196 197 198 200 201 203 204 205 206 357 358 362 364 367 368 370

#cT <- crc[91:100,]# select 50 random bins of 14430
cT <- crc[,3:702]
#cT[,1]="CRC"

#hT <- ctl[91:100,]# select 50 random bins of 23310
hT <- ctl[,3:702]
#hT[,1]="HLT"

cT20 <- cT[,cL]# pick all bins corresponding to the top 20
hT20 <- hT[,cL]# pick all bins corresponding to the top 20

cT1 <- cT[,364]
hT1<-hT[,364]
cT2 <- cT[,205]
hT2<-hT[,205]
cT3 <- cT[,198]
hT3<-hT[,198]

cv364 <- matrix(cT1, ncol = 555, byrow = 26) #
cv205 <- matrix(cT2, ncol = 555, byrow = 26) #  
cv198 <- matrix(cT3, ncol = 555, byrow = 26) #  
indx<-which(kc364>1.5) #  genomic bins: 49  64  90 130 132 209 216 267 308 324 508 519
cv2 <-matrix(,nrow=26,ncol=length(indx))
cv2<-cv364[,indx]
#cv2[,1]<-cv364[,308]
#cv2[,2]<-cv364[,519]
#cv2[,3]<-cv205[,282]
#cv2[,4]<-cv205[,369]
#cv2[,5]<-cv198[,36]
#cv2[,6]<-cv198[,447]
hv364 <- matrix(hT1, ncol = 555, byrow = 42) #  
hv205 <- matrix(hT2, ncol = 555, byrow = 42) # 
hv198 <- matrix(hT3, ncol = 555, byrow = 42) # 
hv2 <-matrix(,nrow=42,ncol=length(indx))
hv2<-hv364[,indx]
#hv2[,1]<-hv364[,308]
#hv2[,2]<-hv364[,519]
#hv2[,3]<-hv205[,282]
#hv2[,4]<-hv205[,369]
#hv2[,5]<-hv198[,36]
#hv2[,6]<-hv198[,447]

dimnames(cv2)[[1]] <- unique(crc[,1]) 
#dimnames(cv2)[[2]] <-  c("308","519","282","369","36","447")

hcTop20 <- rbind(cv2 , hv2)
df<-scale(hcTop20)
col <- colorRampPalette(brewer.pal(11, "RdYlBu"))(256)
hm <- heatmap(df, scale = "none", col =  col) # 3 CTL mixed with CRC: #36, 55, 44, i.e: CTL #10, 29, 18. 
# hm$rowInd CRC 1-26, CTL 27-68 
#(bottom) 25 24 10 22  4 23 44 61 20  6 15 21 11 12 55 13 14 68 36 52  7 19 16  5  8  3  2 17  9 18  1 54 48 51 63 
# 67 28 56 49 62 50 53 33 58 47 45 66 29 30 26 64 46 59 57 42 60 65 27 31 34 32 40 41 38 43 37 35 39 (top end)
# hm$colInd order of 12 bins: 4  6 10  8  7  9  3  2 12  5  1 11

distr <- dist(df)
hcr <- hclust(distr)
hcr$order # 38 41 37 43 35 39 27 31 34 32 40 50 53 47 45 66 33 58 62 49 56 28 67 63 51 48 54 29 30 26 64 57 42 60 65 46 59  
#8  5 16 19 17  2  3  1  9 18 14 13 55 11 12 10 24 25 23 4 22 15 21  6 20 44 61  7 52 36 68

d2 <- dist(hcTop20,method = "euclidean", diag = FALSE, upper = TRUE)
c2 <- hclust(d2, method = "ward.D2", members = NULL)

#heatmap(x, Rowv = NULL, Colv = if(symm)"Rowv" else NULL,
#        distfun = dist, hclustfun = hclust,
#        reorderfun = function(d, w) reorder(d, w),
#        add.expr, symm = FALSE, revC = identical(Colv, "Rowv"),
#        scale = c("row", "column", "none"), na.rm = TRUE,
#        margins = c(5, 5), ColSideColors, RowSideColors,
#        cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc),
#        labRow = NULL, labCol = NULL, main = NULL,
#        xlab = NULL, ylab = NULL,
#        keep.dendro = FALSE, verbose = getOption("verbose"), …)

heatmap.2(df)

y <- df
## Row- and column-wise clustering
hr <-hclust(as.dist(1-cor(t(y), method="pearson")), method="complete")
hc <-hclust(as.dist(1-cor(y, method="spearman")), method="complete")
## Tree cutting
mycl <-cutree(hr, h=max(hr$height)/1.5); mycolhc <-rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)]
## Plot heatmap
mycol <-colorpanel(40, "darkblue", "yellow", "white") 
heatmap.2(y, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=mycol, scale="row", density.info="none")




dNc= c("CRC1","CRC2","CRC3","CRC4","CRC5","CRC6","CRC7","CRC8","CRC9","CRC10")#,"CRC11","CRC12","CRC13","CRC14","CRC15","CRC16",
       #"CRC17","CRC18","CRC19",
      #"CRC20","CRC21","CRC22","CRC23","CRC24","CRC25","CRC26","CRC27","CRC28","CRC29","CRC30","CRC31","CRC32","CRC33","CRC34",
      #"CRC35","CRC36","CRC37","CRC38","CRC39",
      #"CRC40","CRC41","CRC42","CRC43","CRC44","CRC45","CRC46","CRC47","CRC48","CRC49","CRC50","CRC51")

cpts <-unique(crc[,1])
cpts <- unique(crc$sample)[26]

cPtL <- crc [,1] == unique(crc$sample)[1]

crcP1 <- crc[,crc[,1]==cpts[1]]




dimnames(cT20)[[1]] <- unique(crc[,1]) # 27 CRC patients id numbers

c20[3,]<-c20[27,]# 3rd CRC is all zeros
c20<-c20[1:26,]

h20[9,]<-h20[43,]# 9th Healthy is all zeros
h20<-h20[1:42,]

hcTop20 <- rbind(c20 , h20)

#c21<-transpose(c20)# dimension names get confused

df<-scale(hcTop20)
#df<-scale(cM364)
#df<-scale(hM364)


col <- colorRampPalette(brewer.pal(11, "RdYlBu"))(256)

heatmap(df, scale = "none", col =  col)

cT1<-c20[,17]
hT1<-h20[,17]
#rocCH <- append(c20,h20)
TP <- matrix(rep(1, 26*12), ncol = 12, byrow = 26)
TN <- matrix(rep(0, 42*12), ncol = 12, byrow = 42)
rocTF <- rbind(TP , TN)

#rocTF <- append(rep(1, 26),integer(42))
#pred <- prediction(rocCH, rocTF)

pred <- prediction(hcTop20, rocTF)
perf<-performance(pred,"tpr", "fpr")
plot(perf)
auc<- performance(pred,"auc")
auc

sName = bFr[,1] # list of sample names
typeof(sName) # [1] "character"
#list
#typeof(aux)
#[1] "character
#as.numeric(str_extract_all(string, aux)[[1]])
#as.numeric(strsplit(string, "aux")[[1]][-1])

cM364[,308]
#  5  7  8 12  5 12 17  2  7  8 20 17 13 17  7  1  9 12  1 12  8 10  6 10  5 25
hM364[,308]
# 31 20 28 20 33 23 22 19 41 15 33 41 41 28 34 22 35 13 18 22 21 18 14 19 19  5 23 14 15 17 22 22 24 16 11 19 22 23 18 21 13 15
mean(cM364[,308])
# 9.846154
mean(hM364[,308])
# 22.14286
mean(cM364[,519])
# 8.538462
mean(hM364[,519])
# 20.85714

cM364[,519]
# 15 15  8  7  4 12  8  6  8  6  7 11 11  5 13  2  9 11  4  9  6  3 11  8  8 15
hM364[,519]
# 23 16 16 17 31 23 17 42 31 18 43 34 49 34 35 17 29 16 20 14 21 16 22 26 16 10 20 10 10 18  6 17 18 17 15 15 16 13 18 23 10 14

cM198[,447]
# 132 215 186 227 121 294 201 114 226 175 189 193 182 194 163  72 107 192 110 185 197 129 185 210 162 255
hM198[,447]
# 422 365 323 305 426 365 277 431 508 208 412 431 371 379 440 361 535 203 226 251 379 249 274 385 282 214 359 284 200 293 255 226 254 282 145 232 305 216 284 325 246 184
wasserstein1d(cM364[,519], hM364[,519])
wasserstein(cM364,hM364,p=1)