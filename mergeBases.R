library(dplyr)
library(readxl)
library(abind)
#setwd ('~/genomedk/projects/test/KLD_W4')
#Read output from ctDNAtool (Elias) - its my personal folder.
#raw <- read.table("~/genomedk/projects/test/fraglen/out.tsv", sep = "\t", header = F)
PON <- list.files("~/genomedk/projects/test/pon38", recursive = T, full.names = T, pattern = "txt")
PONname <- gsub(".+/(Donor.+)_consensus.+", "\\1", PON)

PON_m <- array(0, c(45,(72376),699))
str(PON_m)

for (i in 1:1){
  i=1
  #raw <- read.table(PON[i], header = F, nrows = 200, stringsAsFactors = F)  
  raw <- read.table(PON[i], header = F, stringsAsFactors = F)  
  
  #res <-
  #  dplyr::group_by(raw, V1, V2) %>%
  #  dplyr::summarise_at(.vars = names(.)[4:ncol(.)], .funs = sum)
  
  auxPON <- raw[,3:701]
  
  PON_m[i,,] <- auxPON
  #write.csv(res,paste0('~/genomedk/projects/test/pon38/',PONname[i],'.csv'))
}
#######################################
PON <- list.files("~/genomedk/projects/test/pon38", recursive = T, full.names = T, pattern = "txt")
raw1 <- read.table(PON[1], header = F, stringsAsFactors = F)  
write.table(A,"/home/matov/projects/test/meanPON38.csv",sep="\t",row.names=TRUE)

meanPON <- read.csv(file = '~/genomedk/projects/test/meanPON38a.csv')
medPON <- meanPON[,2:700]
write.csv(medPON,paste0('~/genomedk/projects/test/pon38/medPON38.csv'))

medPON <- read.csv(file = '~/genomedk/projects/test/pon38/medPON38.csv')
#mdPON <- array(0, c((72376),702))
#mdPON[,1] <- raw1[,1]
#mdPON[,2] <- raw1[,2]
#mdPON[,3:702] <- medPON[,1:700]

aux <- cbind(raw1[,2],medPON)
mdPON <- cbind(raw1[,1],aux)
names(mdPON) <- names(raw1)
 
write.csv(mdPON,paste0('~/genomedk/projects/test/pon38/medPON38v.csv'))

medPONc <-
  dplyr::group_by(mdPON, V1, V2) %>%
  dplyr::summarise_at(.vars = names(.)[4:ncol(.)], .funs = sum)

write.csv(medPONc,paste0('~/genomedk/projects/test/pon38/medPON38kld.csv'))

medPONkld <- read.csv('~/genomedk/projects/test/pon38/medPON38kld.csv')
dim(medPONkld)
######################################
PON_m <- array(0, c(45,(72376/4),699))
str(PON_m)

for (i in 1:2){
  i=1
#raw <- read.table(PON[i], header = F, nrows = 200, stringsAsFactors = F)  
raw <- read.table(PON[i], header = F, stringsAsFactors = F)  

res <-
  dplyr::group_by(raw, V1, V2) %>%
  dplyr::summarise_at(.vars = names(.)[4:ncol(.)], .funs = sum)

auxPON <- res[,3:701]

PON_m[i,,] <- auxPON
#write.csv(res,paste0('~/genomedk/projects/test/pon38/',PONname[i],'.csv'))
}
#mads try:
tmp <- 
  lapply(PON[1:45], function(x){
    as.matrix(read.table(x, header = F, stringsAsFactors = F)[,4:702] )# %>% 
      # dplyr::group_by(V1, V2) %>%
      # dplyr::summarise_at(.vars = names(.)[4:ncol(.)], .funs = sum)
  })
str(tmp)
#Reduce
aa <- simplify2array(tmp) 
A <- apply(aa, c(1,2), median)
head(A)

a1<-matrix(1:4, nrow = 2)
a2<-matrix(10:13, nrow = 2)

A2 <- abind::abind(a1, a2, along = 3) 
apply(A2, c(1,2), sum)

rm(list = ls()) # clear environment

PON_t <- PON_test[,3:701]

PON_test <- read.csv('~/genomedk/projects/test/pon38/PON_merged45.csv')
dim(PON_test)
##########################
QIA <- list.files("~/genomedk/projects/test/qiagen", recursive = T, full.names = T, pattern = "txt")

for (i in 1:(length(QIA))){
  raw2 <- read.table(QIA[i], header = F)  
  
  res2 <-
    dplyr::group_by(raw2, V1, V2) %>%
    dplyr::summarise_at(.vars = names(.)[4:ncol(.)], .funs = sum)
  
  write.csv(res2,paste0('~/genomedk/projects/test/qiagen/QIA_merged',i,'.csv'))
}
QIA_test <- read.csv('~/genomedk/projects/test/qiagen/QIA_merged23.csv')
QIA_t <- QIA_test[,4:702]
write.csv(QIA_t,paste0('~/genomedk/projects/test/KLD_w4/QIAtest.csv'))
dim(QIA_test)
###########################
CRUK <- list.files("~/genomedk/projects/test/fraglen", recursive = T, full.names = T, pattern = "txt")
CRUKname <- gsub(".+/(.+)_consensus.+", "\\1", CRUK)

for (i in 1:(length(CRUKname))){
  raw1 <- read.table(CRUK[i], header = F)  
  
  res1 <-
    dplyr::group_by(raw1, V1, V2) %>%
    dplyr::summarise_at(.vars = names(.)[4:ncol(.)], .funs = sum)
  #res2 <- res1[,4:702]
  write.csv(res1,paste0('~/genomedk/projects/test/fraglen/',CRUKname[i],'.csv'))
}
CRUK_test <- read.csv('~/genomedk/projects/test/fraglen/CRUK_merged123.csv')
dim(CRUK_test)
################################
PON_test <- read.csv('~/genomedk/projects/test/pon38/PON_merged45.csv')
CRUK_test <- read.csv('~/genomedk/projects/test/fraglen/CRUK_merged123.csv')
PON_t <- PON_test[,4:702]
CRUK_t <- CRUK_test[,4:702]
write.csv(PON_t,paste0('~/genomedk/projects/test/KLD_w4/PONtest.csv'))
write.csv(CRUK_t,paste0('~/genomedk/projects/test/KLD_w4/CRUKtest.csv'))
#write.csv(PON_t,'~/genomedk/projects/test/KLD_w4/PONtest.csv')
##############################################
KLDc <- read.csv('~/genomedk/projects/test/KLD_w4/KLD_PONtest_CRUKtest.csv')
plot(KLDc)
KLD <- read.csv('~/genomedk/projects/test/KLD_w4/KLD_PONtest_QIAtest.csv')
plot(KLD)

for (i in 1:224){
  CRUK_aux <- read.csv(paste0('~/genomedk/projects/test/fraglen/CRUK_merged',i,'.csv'))
  CRUK <- CRUK_aux[,4:702]
  write.csv(CRUK,paste0('~/genomedk/projects/test/KLD_w4/CRUK_merged',i,'B.csv'))
}
#x  <- seq(20, 1, 310)
plot(KLDc[20:310,2], type = "l", col = "red")
par(new=TRUE)
plot(KLD[20:310,2], type = "l", col = "green")
#################################################
#KLDc <- read.csv('~/genomedk/projects/test/KLD_w4/KLD_PONtest_CRUK_merged2B.csv')
#KLDc <- read.csv('~/genomedk/projects/test/KLD_w4/KLD_PONtest_CRUKtest.csv')

KLDc <- read.csv('~/genomedk/projects/test/KLD_w4/KLD_PONtest_CRUK_merged114B.csv')
KLD <- read.csv('~/genomedk/projects/test/KLD_w4/KLD_PONtest_QIAtest.csv')
plot(KLDc[20:310,2], ylim=range(c(0,14)), type = "p", col = "red")
par(new=TRUE)
plot(KLD[20:310,2], ylim=range(c(0,14)), type = "p", col = "green")
