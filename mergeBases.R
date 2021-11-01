library(dplyr)
library(readxl)
setwd ('~/genomedk/projects/test/KLD_W4')

#Read output from ctDNAtool (Elias) - its my personal folder.
#raw <- read.table("~/genomedk/projects/test/fraglen/out.tsv", sep = "\t", header = F)
PON <- list.files("~/genomedk/projects/test/pon38", recursive = T, full.names = T, pattern = "txt")

PONname <- gsub(".+/(Donor.+)_consensus.+", "\\1", PON)

for (i in 1:45){
  i=1
#raw <- read.table(PON[i], header = F, nrows = 200, stringsAsFactors = F)  
raw <- read.table(PON[i], header = F)  
res <-
  dplyr::group_by(raw, V1, V2) %>%
  dplyr::summarise_at(.vars = names(.)[4:ncol(.)], .funs = sum)

#write.csv(res,paste0('~/genomedk/projects/test/pon38/',PONname[i],'.csv'))
}

PON_t <- PON_test[,4:702]

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

#########################

###########################
CRUK <- list.files("~/genomedk/projects/test/fraglen", recursive = T, full.names = T, pattern = "txt")


for (i in 1:(length(CRUK))){
  raw1 <- read.table(CRUK[i], header = F)  
  
  res1 <-
    dplyr::group_by(raw1, V1, V2) %>%
    dplyr::summarise_at(.vars = names(.)[4:ncol(.)], .funs = sum)
  
  write.csv(res1,paste0('~/genomedk/projects/test/fraglen/CRUK_merged',i,'.csv'))
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

#x  <- seq(20, 1, 310)
plot(KLDc[20:310,2], type = "l", col = "red")
par(new=TRUE)
plot(KLD[20:310,2], type = "l", col = "green")
