library("dplyr")
library("glmnet")
source("~/genomedk/matovanalysis/umiseq_analysis/R/read_bed.R") #1/0 list
# add new flags; number of fragments, age, gender, concentration, 10 flags
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
v0 <- min(v[v>0])/10000000 # for counts w zero variance, we replace w a very small value
v1<- v
v1[v==0]=v0
# PON mutations and variability ###################################################################################
sitemut <- t(apply(pon_obj2$coordinates, 1, function(x){
  paste( paste0(trimws(x[1]), ":", trimws(x[2]), "_"),
         paste(x[3] ,c("A", "T", "C", "G"), sep = "/"),
         sep = "")}))
#dfS <- data.frame(sitemut) 
#sitemutPON <- dfS[list == 1,] 
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
    mahaT[j] <- sum(countsW[i,,]* mafsC[k,,]/v1,  na.rm=T)
    j = j + 1
  }
  print(mahaT)
}
length(mahaT)
plot(mahaT)

mahaT2 <- vector()
j=1
for (i in cancer_list){
  #i = 6
  for (k in adenoma_list){
    mahaT2[j] <- sum(countsW[i,,]* mafsC[k,,]/v1,  na.rm=T)
    j = j + 1
  }
  print(mahaT2)
}
length(mahaT2)
plot(mahaT2)

maha <- vector()
j=1
for (i in cancer_list){
maha[j] <- sum(countsW[i,,]* mafsC[i,,]/v1,  na.rm=T) # should be VAF mafsC[i,,]) #45 of 95 cancers have WES, of them 38 have VAF>0
# maha[i] <- sum(countsW[6,,]* mafsC[i,,]/v1,  na.rm=T) # should be VAF mafsC[i,,]) #
  if (maha[j]>0) {
j = j + 1
}
#aux <- sum(countsW[6,,]*countsC001[i,,]/v1,  na.rm=T) # should be VAF mafsC[i,,]) #
# for countsC001 #6 [1] 142
print(maha)
# devide by PON var each pos
}
plot(maha) # Cancer samples with VAF>0
mahaCancer <- maha[1:37] # 37 of 45 , w 8 cancers w VAFs 0
mahaControl <- c(mahaT, mahaT2) # 1575 (45x15+45x20)

#ROC; find where overall success rate numbers are
#condition <- rbind(array(1, dim=c(79,12)), array(0, dim=c(74,12)))
#pred <- prediction(c(mahaCancer, mahaControl), condition, label.ordering = c(0, 1))  
rocTF = rep(0, (37+45*15+45*20))
rocTF[1:37]=1
pred <- prediction(c(mahaCancer, mahaControl), rocTF)
perf<-performance(pred,"tpr", "fpr")
plot(perf)
auc<- performance(pred,"auc")
auc
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


#maha1 <- vector()
#maha1 [1:69]<- maha[2:70] 
#maha1 [70:95] <- maha[105:130]
#plot(maha1)


######################################################################################
mahaAll <- vector()
for (i in 1:130){
  maha[i] <- sum(countsW[i,,]* mafsC[i,,]/v1,  na.rm=T) # should be VAF mafsC[i,,]) #
  
  #aux <- sum(countsW[6,,]*countsC001[i,,]/v1,  na.rm=T) # should be VAF mafsC[i,,]) #
  # for countsC001 #6 [1] 142
  print(maha[i])
  # devide by PON var each pos
}
plot(mahaAll)
mahaAllCancer <- vector()
mahaAllCancer [1:69]<- mahaAll[2:70] 
mahaAllCancer [70:95] <- mahaAll[105:130]
plot(mahaAllCancer)

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

# extract the sample name from the first pileup list in dimnames


j=1
k = 3
for (i in 1:k){
i=1
aux <- mafsC1[j,,] # plasma sample 
which(grepl(plasma$library_id[r], pileups)), , ]
mafScore <- aux[sitemutPON == chr5:112815487]  #cancer_SNPs$sitemut_hg38[i]]
scores
}
#####################################################################################
countsW <- readRDS("~/genomedk/matovanalysis/umiseq_analysis/R/cruki.RDS") # 




# IMPROVE WES data ###################################################################################
wes <- read.table("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/umiseq_paper/data/201123_wes-spora-mutations-improve.csv", header = TRUE) # HG38
pts <- unique(wes$pt_id)
# IMPROVE clinical data ###################################################################################
ItW <- read.table("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/umiseq_paper/data/IMPROVEptList",header = TRUE)
ItW$index <- sapply(as.character(ItW$library_id), function(x) grep(x, pileupsIw))  
#################################################################
wes_id <- ItW$pt_id %in% pts
preop_i <- ItW[(ItW$op_time_cat == -1)&wes_id, "index"] # 57 of which 44 have WES
posop14_i <- ItW[ (ItW$op_time_cat == 2)&wes_id, "index"] # 42
posop30_i <- ItW[ (ItW$op_time_cat == 30)&wes_id, "index"] # 42
countsIw <-  piles_to_counts(files = pileupsIw[preop_i], regions = pon_obj2$regions) # PREOP, POSTOP14, POSTOP30
countsW <- countsIw[,,1:4] + countsIw[,,6:9]
mafsW= array(0, dim=c(dim(countsW)[1],dim(countsW)[2],dim(countsW)[3]))
auxM <- rowSums(countsW, dims = 2) 
for (i in 1:dim(countsW)[1]) {
  i=1
  mafsW[i,,] <- countsW[i,,]  /auxM[i,]
}
#################################################################
ww = array(0, dim=c(dim(co)[1],dim(co)[2],dim(co)[3]))
sc<- vector()
for (i in 1:dim(mafsW)[1]) { 
  for (j in 1:length(pts)) {  
    if (grepl(as.character(pts[j]), as.character(pileupsIw[preop_i][i]))) { # for pts with WES data
      wesP <- wes[wes==pts[j],] # find the mutations in WES
      mu <- wesP$sitemut_hg38 # look up the indexes in the sitemut panel
      for (k in 1:length(mu)) { # computer for all WES mutations for this sample
        auxW <- mafsW[i,,]
        iW <- sitemut == mu[k] 
        ww[k] <- auxW[iW]/vo1[iW]*prior11[iW] # VAF divided by PON variance & weighted by Cosmic
      }
      sc[j] <- sum(ww)/length(mu) # mutation score per patient
    }
  }
}
###############plot mutation score########################################################
r1<- min(log2(sw))-1
r2 <- max(log2(sw))+1
plot(log2(sw), ylim=range(c(r1,r2)), col="green", pch = 17)
###### debug ########
#sapply(as.character(pts), function(x) grep(x, pileupsIw[preop_i])) 
#wes_id <- intersect(ItW$pt_id, pts)
#ItW$index <- sapply(as.character(ItW$library_id), function(x) grep(x, pileupsIw)) 
########## mutations for PON1 subject ################
wesP <- wes[wes==pts[1],]
m1 <- wesP$sitemut_hg38[1]
#chr17:7675235_T/C
#chr3:179218294_G/A
#chr7:140753336_A/T
sum(sitemut == "chr17:7675235_T/C") # 1
sum(sitemut == m1) # 1
which(sitemut == m1)# 50720
sitemut[which(sitemut == m1)] # "chr17:7675235_T/C"

########VERSIO 2########################################################################
### Intersect WES mutations with plasma and calculate score ###-----------------
library("dplyr")

wes <- read.table("~/genomedk/N140_Targeting/specs/umiseq_paper/data/201123_wes-spora-mutations-improve.csv", header = TRUE) # HG38
plasma <- read.table("~/genomedk/N140_Targeting/specs/umiseq_paper/data/IMPROVEptList", header = TRUE)
prior <- readRDS("~/genomedk/N140_Targeting/specs/specs_analysis/sw_input_files/180903_prior.RDS")[, 1:4]
pon38 <- readRDS("~/genomedk/N140_Targeting/specs/umiseq_paper/reference/201020_hg38-novaseq-xgen-sporacrc-pon.RDS")
pon19 <- readRDS("~/genomedk/N140_Targeting/specs/umiseq_paper/reference/200419_novaseq-xgen-sporacrc-pon.RDS")

#Add sitemut_hg38 to plasma table (overwrite plasma)
plasma <- dplyr::left_join(plasma, wes[, c("pt_id", "sitemut_hg38")])
plasma <- plasma[!is.na(plasma$sitemut_hg38), ] #Exclude samples with no sitemut


#Add the index for the 18094*4 panel matrix
sitemut <- t(apply(pon38$coordinates, 1, function(x){
  paste( paste0(trimws(x[1]), ":", trimws(x[2]), "_"),
         paste(x[3] ,c("A", "T", "C", "G"), sep = "/"),
         sep = "")}))

plasma$i <- sapply(as.character(plasma$sitemut_hg38), function(s) which(sitemut %in% s))

#Get the pileups
pileups <- list.files("~/genomedk/IMPROVE/sporacrc/N227", recursive = T, full.names = T, pattern = "consensus.bait.pileup$")
pileups <- unique(pileups[sapply(plasma$library_id, function(l)which(grepl(l, pileups)))]) #Only get pileups needed

#Make counts assuming hg19 (otherwise, split for hg19 and hg38)
tmp <- piles_to_counts(pileups, pon19$regions)
counts <- tmp[,,1:4] + tmp[,,6:9]
mafs <- abind::abind( lapply(1:dim(counts)[[1]], function(s)(counts[s,,]/rowSums(counts[s,,]))), along = 0)
vars <- apply(counts, 2:3, var) # variance of the counts of each position based on PON
vars[ vars == 0 ] <- 1e-6


#Calculate f1 per sitemut (row in plasma)
f1 <- function(maf, pon_var, prior)prior * maf/pon_var # VAF divided by PON variance & weighted by Cosmic

plasma$score <-
  sapply(1:nrow(plasma), function(r) {
    maf0 <- mafs[ which(grepl(plasma$library_id[r], pileups)), , ][ plasma$i[r] ]  
    pon_var0 <- vars[ plasma$i[r] ]
    prior0 <- prior[ plasma$i[r] ]
    f1(maf0, pon_var0, prior0 )
  })

#Then sum score for each library_id (overwrite plasma)
plasma <-
  dplyr::group_by(plasma, library_id) %>%
  dplyr::mutate(sum = sum(score))