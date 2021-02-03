source("~/genomedk/matovanalysis/umiseq_analysis/R/read_bed.R") #1/0 list
setwd ('~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/specs_analysis')
source("sw_input_files/duplex_tools.R")
library("dplyr")
library("ggpubr")
pon_hg19 <- readRDS("~/genomedk/PolyA/faststorage/BACKUP/IMPROVE/call/references/200419_novaseq-xgen-sporacrc-pon.RDS") # 

f3 <- function(a, M, S){
  which(apply(a, c(2, 3), function(x)sum(sum(x >= M, na.rm=T) >= S))>0)
}
S =5; M=0.01;  
#indP = f3(mafsP2, M, S)
indP = c(1598,2207,5434,8551,9245,11243,12551,13815,18170,18937,21378,21687,28047,29455,32718,35119,35254,37031,37182,39129,41149,47691,47769,49058,49983,50282,50377,51134,51276,51280,51289,54657,54661,54990,55624,63322,63689,64088,64382,66794,66987,68940)
#indP = c(1029,1638,4717, 7834, 9768, 11594, 18341, 27868, 33302, 35174, 40763, 42504, 42803, 42898, 43655, 43797, 43801, 43810, 46627, 46631, 47168, 54718, 58832)
VAFcut = 0.35
# 45 Subjects of the Control Panel of Normal PON ####################################################################
pon_obj2 <- readRDS("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/umiseq_paper/reference/201217_hg38-novaseq-xgen-sporacrc-pon.RDS") # 46
# pon_obj2 <- readRDS("sw_input_files/201020_hg38-novaseq-xgen-sporacrc-pon.RDS") # 45 subjects only, 1 is missing
pon_counts <- pon_obj2[["pon"]]

no0 = array(0, dim=c(dim(pon_counts)[1]-1,dim(pon_counts)[2],dim(pon_counts)[3]))
no0[1:27,,] <- pon_counts[1:27,,]
no0[28:45,,]<-pon_counts[29:46,,]
#no0<-pon_counts

no1 = array(0, dim=c(dim(no0)[1],sum(list),dim(no0)[3]))
for (i in 1:dim(no0)[1]) {
  #i=1
  no0[i,,][indP]<-NA # blacklist
  p2 <- data.frame(no0[i,,])#PON
  p1 <- p2 [list == 1, ] 
  no1[i,,] <- data.matrix(p1)
  #no1[i,,][indP]<-NA # blacklist
}
no <- no1[,,1:4]+no1[,,6:9]
mafsP1 = array(0, dim=c(dim(no)[1],dim(no)[2],dim(no)[3]))
erP1<- vector()
auxMP <- rowSums(no, dims = 2) 
for (i in 1:dim(no)[1]) {
  #i = 1
  mafsP1[i,,] <- no[i,,]  /auxMP[i,]
  erP1[i] <- mean(mafsP1[i,,][mafsP1[i,,]<= VAFcut], na.rm=T)  
}
plot(erP1)
print(erP1)

#oL <- erP1[28]
#erP<-vector()
#erP[1:27]<- erP1[1:27]
#erP[28:45] <- erP1[29:46]
#mean(erP) + 3*sd(erP)
#mean(erP1) + 3*sd(erP1)
#oL 
#shapiro.test(erP1)
#shapiro.test(erP)

#boxplot(erP1)
#boxplot(erP)
#boxplot(erQ1)

er <- list(errPno28=erP1, errQ17=erQ1, errI229=erI1, errC69=erC1[CRUKages>1], errDS8=erD1)
boxplot(er,notch = TRUE,horizontal = TRUE,border = "brown",col = c("green","blue","red","orange","black"))

PON<- read_xlsx('~/genomedk/matovanalysis/umiseq_analysis/2020-11-04_PON_age_gender.xlsx')
PON1 <- list()
PON1<-PON[1:27,]
PON11<-rbind(PON1,PON[29:46,])
age <- list(PON_age =PON11$age, F_age = PON11$age[PON11$gender=="F"], M_age = PON11$age[PON11$gender=="M"] )
boxplot(age,notch = TRUE,horizontal = TRUE,border = "brown",col = c("blue", "green", "red"))

er1 <- list(ePover53=erP1[PON11$age>53], ePless54=erP1[PON11$age<54])
boxplot(er1,notch = TRUE,horizontal = TRUE,border = "brown",col = c("blue","red"))

PON1 <- list()
PON1<-PON[1:27,]
PON11<-rbind(PON1,PON[29:46,])

p <- ggboxplot(ErrorRates, x = "gender", y = "age",
               color = "gender", palette = "jco",
               add = "jitter")
p + stat_compare_means()
p + stat_compare_means(method = "t.test")
IMPROVE <- read_xlsx('~/genomedk/matovanalysis/umiseq_analysis/2020-01-21_IMPROVE_samples_age.xlsx')

erI1 #229
PON11$age>53


CRUK <- read_xlsx('~/genomedk/matovanalysis/umiseq_analysis/2020-01-04_CRUK_sample_status.xlsx')
CRUKlist <- CRUK$`Biobank label`[CRUK$sequenced=="yes"]#80

listCRUK <- unlist(sapply(CRUKlist, function(x) grep(x, x = pileupsC )))
pileupsC[listCRUK] # 69

CRUK$pt_age[listCRUK]
CRUKages <- CRUK$pt_age[listCRUK] # 69
CRUKages>70 & CRUKages<81
erC1[CRUKages>80] # 10
erC1[CRUKages>70 & CRUKages<81] # 25
erC1[CRUKages>60 & CRUKages<71] # 19
erC1[CRUKages>50 & CRUKages<61] # 11
erC1[CRUKages<51] # 4


length(IMPROVE$age_at_sample_time)
length(erI1)

erCI <- list(Int5=c(erC1[CRUKages>80],erI1[IMPROVE$age_at_sample_time>80]), 
             Int4=c(erC1[CRUKages>70 & CRUKages<81],erI1[IMPROVE$age_at_sample_time>70 & IMPROVE$age_at_sample_time<81]),
             Int3=c(erC1[CRUKages>60 & CRUKages<71],erI1[IMPROVE$age_at_sample_time>60 & IMPROVE$age_at_sample_time<71]),
             Int2=c(erC1[CRUKages>50 & CRUKages<61],erI1[IMPROVE$age_at_sample_time>50 & IMPROVE$age_at_sample_time<61]),
             Int1=c(erC1[CRUKages<51],erI1[IMPROVE$age_at_sample_time<51]))
boxplot(erCI,notch = TRUE,horizontal = TRUE,border = "brown",col = c("red","orange", "brown","blue","green"))

listCRUKhe <- unlist(sapply(CRUK$`Biobank label`, function(x) grep(x, x = pileupsC[1:8])))#7
CRUK$pt_age[CRUK$`Biobank label`=="S07A05776D"]#"53"
erCH <- 2.173574e-05  
CRUKheAge <- 53

erP1[PON11$age<54]

QIAGEN<- read_xlsx('~/genomedk/matovanalysis/umiseq_analysis/2020-11-04_PON_age_gender.xlsx', sheet = 2)
qiagenList <- list()
qiagenList[1:15]<-pileupsQ[1:15]
qiagenList[16:22]<-pileupsQ[18:24]
listQiagen <- unlist(sapply(QIAGEN$donor, function(x) grep(x, x = unlist(qiagenList))))#7

#16 and 17 are D415 (#22 in the list) and D1416 (#7 in the list)
QIAGEN$age
qAge <- vector()
qAge[1:15]<-QIAGEN$age[1:15]
qAge[16:22]<-QIAGEN$age[18:24]

erPQ <- list(Int3=c(erQ1[qAge>60 & qAge<71],erP1[PON11$age>60 & PON11$age<71]),
             Int2=c(erQ1[qAge>50 & qAge<61],erP1[PON11$age>50 & PON11$age<61],erCH),
             Int1=c(erQ1[qAge<51],erP1[PON11$age<51]))
boxplot(erPQ,notch = TRUE,horizontal = TRUE,border = "brown",col = c("brown","blue","green"))

# QIAGEN healthy samples ############################################################################################
pileupsQ <- list.files("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/qiagen_kit_test/201019", recursive = T, full.names = T, pattern = "bait.pileup")
countsQ00 <-  piles_to_counts(files = pileupsQ, 
                              regions = pon_hg19$regions)
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


### IMPROVE #####
pileupsI <- list.files("~/genomedk/PolyA/faststorage/BACKUP/IMPROVE/sporacrc/N227", recursive = T, full.names = T, pattern = "bait.pileup$")
countsI0 <-  piles_to_counts(files = pileupsI[1:213], # DEC 8 noon, there are 61 IMPROVE files.
                             regions = pon_hg19$regions)
countsI= array(0, dim=c(dim(countsI0)[1],sum(list),dim(countsI0)[3]))
for (i in 1:dim(countsI0)[1]) {
  countsI0[i,,][indP]<-NA # blacklist
  p2 <- data.frame(countsI0[i,,]) 
  p1 <- p2 [list == 1, ] 
  countsI[i,,] <- data.matrix(p1)
}
countsI1 <- countsI[,,1:4] + countsI[,,6:9]
erI1<- vector()
mafsI1 = array(0, dim=c(dim(countsI1)[1],dim(countsI1)[2],dim(countsI1)[3]))
auxMI <- rowSums(countsI1, dims = 2) 
for (i in 1:dim(countsI1)[1]) {
  mafsI1[i,,] <- countsI1[i,,]  /auxMI[i,]
  erI1[i] <- mean(mafsI1[i,,][mafsI1[i,,]<= VAFcut], na.rm=T) 
}
plot(erI1)
# CRUK patient samples ##########################################################################################
pileupsC <- list.files("~/genomedk/PolyA/faststorage/BACKUP/CRUK/plasma/N289", recursive = T, full.names = T, pattern = "bait.pileup")
countsC0 <-  piles_to_counts(files = pileupsC[1:90], #pileupsC[listCRUK], #, 
                             regions = pon_hg19$regions)
countsC= array(0, dim=c(dim(countsC0)[1],sum(list),dim(countsC0)[3]))
for (i in 1:dim(countsC0)[1]) {
  countsC0[i,,][indP]<-NA # blacklist
  p2 <- data.frame(countsC0[i,,])#CRUK
  p1 <- p2 [list == 1, ] 
  countsC[i,,] <- data.matrix(p1)
}
countsC1 <- countsC[,,1:4] + countsC[,,6:9]
erCH1<- vector()
mafsC1 = array(0, dim=c(dim(countsC1)[1],dim(countsC1)[2],dim(countsC1)[3]))
auxMC <- rowSums(countsC1, dims = 2) 
for (i in 1:dim(countsC1)[1]) {
  mafsC1[i,,] <- countsC1[i,,]  /auxMC[i,]
  erCH1[i] <- mean(mafsC1[i,,][mafsC1[i,,]<= VAFcut], na.rm=T) 
}
plot(erCH1)
# DS samples ##################################################################################################
dat0 <- readRDS("sw_output_files/2020-10-23-145546_sw-output.RDS")  
pileupsD <- unlist(attributes(dat0))
pileupsD <- sub("/faststorage/project/PolyA/BACKUP", "~/genomedk/PolyA/faststorage/BACKUP", pileupsD)
all(file.exists(pileupsD))
countsD <-  piles_to_counts(files = pileupsD, 
                            regions = pon_obj2$regions)
countsDD= array(0, dim=c(dim(countsD)[1],sum(list),dim(countsD)[3]))

founderinfo <- readRDS("200330_founder-mutations-method1.RDS")
fI<-founderinfo[1:90,]# only the 90 mutations info
inMu <- fI$index# indexes of the 90 mutations on the panel

for (i in 1:dim(countsD)[1]) {
  countsD[i,,][inMu]<-NA # panel
  countsD[i,,][indP]<-NA # blacklist
  p2 <- data.frame(countsD[i,,])
  p1 <- p2 [list == 1, ] 
  countsDD[i,,] <- data.matrix(p1)
}
countsD1 <- countsDD[,,1:4] + countsDD[,,6:9]
mafsD1 = array(0, dim=c(dim(countsD1)[1],dim(countsD1)[2],dim(countsD1)[3]))
erD1<- vector()
auxMD <- rowSums(countsD1, dims = 2) 
for (i in 1:dim(countsD1)[1]) {
  mafsD1[i,,] <- countsD1[i,,]  /auxMD[i,]
  erD1[i] <- mean(mafsD1[i,,][mafsD1[i,,]<= VAFcut], na.rm=T)
}
plot(erD1)
###########################################################################################################
plot(erQ1, na.rm=T,ylim=range(c(1e-05,0.00011)), col="green", pch = 17)
par(new = TRUE)
plot(erP1, na.rm=T,ylim=range(c(1e-05,0.00011)), col="purple", pch = 12)
legend(3,1.1e-04,legend=c("QIAGEN", "PON"),col=c("green","blue"),lty=1:1, cex=1.0)
mean(erQ1, na.rm=T)
mean(erP1, na.rm=T)
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