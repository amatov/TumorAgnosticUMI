library("glmnet")
source("~/genomedk/matovanalysis/umiseq_analysis/R/read_bed.R") #1/0 list


ctl1D23 = array(0, dim=c(74,499))
for (i in 1:74  ) {
  print(i)
  #i=2
  auxFR <- read.table(pileupsD2[listCTL1[i]], header = TRUE) # sample per sample, file per file. 
  ctl1D23[i,] <- unlist(colSums(auxFR[,2:500])) / sum(unlist(auxFR[,2:500]))
  
}
colD23 = array(0, dim=c(79,499))
for (i in 1:79  ) {
  print(i)
  #i=2
  auxFR <- read.table(pileupsD2[listCOL[i]], header = TRUE) # sample per sample, file per file. 
  colD23[i,] <- unlist(colSums(auxFR[,2:500]))  / sum(unlist(auxFR[,2:500]))
}
#x <- rbind(ctl1D23[1:37,], colD23[1:37,])
#y = rep(0, 153)
#y[75:153]=1

samplesTr = array(0, dim=c((37+40),499))
samplesTr <- rbind(ctl1D23[1:37,], colD23[1:40,])
selectionTr = rep(0, 77)
selectionTr[38:77]=1

samplesTe = array(0, dim=c((37+39),499))
samplesTe <- rbind(ctl1D23[38:74,], colD23[41:79,])
selectionTe= rep(0, 76)
selectionTe[38:76]=1 

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
################################################################
colD2N = array(0, dim=c(129,574,499))
j=1
for (i in 1:129  ) {
  print(i)
  #i=2
  if (i<80){
  auxFR <- read.table(pileupsD2[listCOL[i]], header = TRUE) # sample per sample, file per file. 
  } else {
    k <- i-79
  auxFR <- read.table(pileupsD2[listREC[k]], header = TRUE) # sample per sample, file per file. 
  }
  colD2N[j,,] <- unlist(auxFR[,2:500])/ sum(unlist(auxFR[,2:500]))
  j=j+1
}
col4D2N = array(0, dim=c(23,574,499))
j=1
for (i in 1:23 ) {
  print(i)
  #i=2
  auxFR <- read.table(pileupsD2[listCOL4[i]], header = TRUE) # sample per sample, file per file. 
  col4D2N[j,,] <- unlist(auxFR[,2:500])/ sum(unlist(auxFR[,2:500]))
  j=j+1
}
# all stages 79 colon cancers vs 74 control1 no comorbidity
samplesTr = array(0, dim=c((37+40),length(ind195)))
#samplesTr <- rbind(ctl1D23[1:37,], colD23[1:40,])
samplesTr <- rbind(ctl1D2N[,ind195,195][1:37,], colD2N[,ind195,195][1:40,])
selectionTr = rep(0, 77)
selectionTr[38:77]=1
samplesTe = array(0, dim=c((37+39),length(ind195)))
#sampplesTe<- rbind(ctl1D23[38:74,], colD23[41:79,])
samplesTe <- rbind(ctl1D2N[,ind195,195][38:74,], colD2N[,ind195,195][41:79,])
selectionTe= rep(0, 76)
selectionTe[38:76]=1 

# only stage IV 23 colon cancers vs 74 control1 no comorbidity
k=1
samplesTr = array(0, dim=c((37+12),(length(ind195)+k)))
samplesTr[,1:length(ind195)] <- rbind(ctl1D2N[,ind195,195][1:37,], col4D2N[,ind195,195][1:12,])
samplesTr[,(length(ind195)+k)] <- c(m2$age[listCTL1][1:37],m2$age[listCOL4][1:12])

# add new flags; number of fragments, age, gender, concentration, 10 flags
selectionTr = rep(0, 49)
selectionTr[38:49]=1
samplesTe= array(0, dim=c((37+11),(length(ind195)+k)))
samplesTe[,1:length(ind195)] <- rbind(ctl1D2N[,ind195,195][38:74,], col4D2N[,ind195,195][13:23,])
samplesTe[,(length(ind195)+k)] <- c(m2$age[listCTL1][38:74],m2$age[listCOL4][13:23])
selectionTe= rep(0, 48)
selectionTe[38:48]=1 

# umiseq 56 PreOp IMPROVE vs 45 PON
samplesTr = array(0, dim=c((23+28),595))
samplesTr <- rbind(umiN[,,356][1:23,], umiiN[,,356][1:28,])
selectionTr = rep(0, 51)
selectionTr[24:51]=1
samplesTe = array(0, dim=c((22+28),595))
samplesTe <- rbind(umiN[,,356][24:45,], umiiN[,,356][29:56,])
selectionTe= rep(0, 50)
selectionTe[23:50]=1 

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

# Splitting the NORMALIZED data into test and train
samplesTr = array(0, dim=c((37+40),574*2))
ctlTr <- rbind(t(ctl1D2N[,,195][1:37,]), t(ctl1D2N[,,137][1:37,]))
#ctlTr <- rbind(ctlTr, t(ctl1D2[,,365][1:37,]))
colTr <- rbind(t(colD2N[,,195][1:40,]), t(colD2N[,,137][1:40,]))
#colTr <- rbind(colTr, t(colD2[,,365][1:40,]))
samplesTr <- rbind(t(ctlTr),t(colTr))# dim(samplesTr) 77 574 w 37 ctl and 40 cc
selectionTr = rep(0, 77)
selectionTr[38:77]=1

samplesTe = array(0, dim=c((37+39),574*2))
ctlTe <- rbind(t(ctl1D2N[,,195][38:74,]), t(ctl1D2N[,,137][38:74,]))
#ctlTe <- rbind(ctlTe, t(ctl1D2[,,365][38:74,]))
colTe <- rbind(t(colD2N[,,195][41:79,]), t(colD2N[,,137][41:79,]))
#colTe <- rbind(colTe, t(colD2[,,365][41:79,]))
samplesTe <- rbind(t(ctlTe), t(colTe))
selectionTe= rep(0, 76)
selectionTe[38:76]=1 

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
