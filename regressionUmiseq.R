library("glmnet")
source("~/genomedk/matovanalysis/umiseq_analysis/R/read_bed.R") #1/0 list

# add new flags; number of fragments, age, gender, concentration, 10 flags

CRUK <- read.table('~/genomedk/matovanalysis/umiseq_analysis/R/specs_data_cruk-plasma-info.lst', header = T)
CRUKlist <- which(CRUK$sample_type=="CRC pre-OP" & CRUK$cancer==1)  
length(CRUKlist) #183 not 130
############################# 
v<-apply(mafsP1,2:3,var) # variance of the counts of each position based on PON
#v<-apply(mafsP1,2:3,var) # variance of the VAFs of each position based on PON
v0 <- min(v[v>0])/10000000 # for counts w zero variance, we replace w a very small value
v1<- v
v1[v==0]=v0
#######################################################################################################################################
#pileupsC <- list.files("~/genomedk/PolyA/faststorage/BACKUP/CRUK/plasma/N289", recursive = T, full.names = T, pattern = "bait.pileup")
#CRUK <- read_xlsx('~/genomedk/matovanalysis/umiseq_analysis/2020-01-04_CRUK_sample_status.xlsx')
#CRUKlist <- CRUK$`Biobank label`[CRUK$sequenced=="yes"]#80

#listCRUK <- unlist(sapply(CRUKlist, function(x) grep(x, x = pileupsC[45:152] )))
#pileupsC[45:152] [listCRUK]  

#countsC0 <-  piles_to_counts(files = pileupsC[45:152] [listCRUK]  , regions = pon_obj2$regions)
countsC00 <- readRDS("~/genomedk/matovanalysis/umiseq_analysis/R/cruk-counts.RDS") # 
countsC0= array(0, dim=c(95,dim(countsC00)[2],dim(countsC00)[3]))
countsC0[1:69,,]<-countsC00[2:70,,]
countsC0[70:95,,]<-countsC00[105:130,,]
normalC0= array(0, dim=c(35,dim(countsC00)[2],dim(countsC00)[3]))
normalC0[1,,]<-countsC00[1,,]
normalC0[2:35,,]<-countsC00[71:104,,]

countsC= array(0, dim=c(dim(countsC0)[1],sum(list),dim(countsC0)[3]))
for (i in 1:dim(countsC0)[1]) {
  #i=1
  p2 <- data.frame(countsC0[i,,])#CRUK
  p1 <- p2 [list == 1, ] 
  countsC[i,,] <- data.matrix(p1)
}
countsC1 <- countsC[,,1:4] + countsC[,,6:9]
mafsC1 = array(0, dim=c(dim(countsC1)[1],sum(list),dim(countsC1)[3]))
mafsC2= array(0, dim=c(dim(countsC1)[1],sum(list)*dim(countsC1)[3]))
auxMC <- rowSums(countsC1, dims = 2) 
for (i in 1:dim(countsC1)[1]) {
  mafsC1[i,,] <- countsC1[i,,]  /auxMC[i,]/v1/sum(list)
  mafsC2[i,] <- mafsC1[i,,]/sum(mafsC1[i,,])#
}
# 45 Subjects of the Control Panel of Normal PON ####################################################################
#pon_obj2 <- readRDS("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/umiseq_paper/reference/201217_hg38-novaseq-xgen-sporacrc-pon.RDS") # 46
#pon_counts <- pon_obj2[["pon"]]
#no0 = array(0, dim=c(dim(pon_counts)[1]-1,dim(pon_counts)[2],dim(pon_counts)[3]))
#no0[1:27,,] <- pon_counts[1:27,,]
#no0[28:45,,]<-pon_counts[29:46,,]
no0 <- normalC0
no1 = array(0, dim=c(dim(no0)[1],sum(list),dim(no0)[3]))
for (i in 1:dim(no0)[1]) {
  #i=1
  p2 <- data.frame(no0[i,,])#PON
  p1 <- p2 [list == 1, ] 
  no1[i,,] <- data.matrix(p1)
}
no <- no1[,,1:4]+no1[,,6:9]
mafsP1 = array(0, dim=c(dim(no)[1],sum(list),dim(no)[3]))
mafsP2 = array(0, dim=c(dim(no)[1],sum(list)*dim(no)[3]))
auxMP <- rowSums(no, dims = 2) 
no2= array(0, dim=c(dim(no)[1],sum(list)*dim(no)[3]))
for (i in 1:dim(no)[1]) {
  mafsP1[i,,] <- no[i,,]  /auxMP[i,]/v1/sum(list)
  mafsP2[i,] <- mafsP1[i,,]/sum(mafsP1[i,,])
}

# umiseq 62k COUNTS 95 PreOp CRUK vs 45 PON
samplesTr = array(0, dim=c((18+48),61860))
samplesTr <- rbind(mafsP2[1:18,], mafsC2[1:48,])
selectionTr = rep(0, 66)
selectionTr[19:66]=1
samplesTe = array(0, dim=c((17+47),61860))
samplesTe <- rbind(mafsP2[19:35,], mafsC2[49:95,])
selectionTe= rep(0, 64)
selectionTe[18:64]=1 

# run #####
cvm = cv.glmnet(samplesTr, selectionTr, family = "binomial", alpha=1, nfolds=10) 

# identifying best lamda
best_lam <- cvm$lambda.min
best_lam 
# Rebuilding the model with best lamda value identified
lasso_best <- glmnet(samplesTr, selectionTr, alpha = 1, lambda = best_lam)
pred <- predict(lasso_best, s = best_lam, newx = samplesTe)

plot(cvm)
plot(cvm$glmnet.fit)

final <- cbind(selectionTe, pred)
final
plot(final)
########################################################################################
plot(final[1:22,2],ylim=c(0.0066,0.0073))
 