library(FactoMineR)
library("RColorBrewer")
library("dplyr")
library(pvclust)
#source("~/genomedk/matovanalysis/umiseq_analysis/R/read_bed.R")
source("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/umiseq_paper/R/read_bed.R")
setwd ('~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/specs_analysis')
source("sw_input_files/duplex_tools.R")

# 45 Subjects of the Control Panel of Normal PON ####################################################################
#pon_obj2 <- readRDS("sw_input_files/201020_hg38-novaseq-xgen-sporacrc-pon.RDS") # 
pon_obj2 <- readRDS("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/umiseq_paper/reference/201217_hg38-novaseq-xgen-sporacrc-pon.RDS") # 46
no1 <- pon_obj2[["pon"]] # counts
no <- no1[,,1:4] + no1[,,6:9]
no2 = array(0, dim=c(dim(no)[1],dim(no)[2]*dim(no)[3]))
mafsP1 = array(0, dim=c(dim(no)[1],dim(no)[2],dim(no)[3]))
auxMP <- rowSums(no, dims = 2) 
for (i in 1:dim(no)[1]) {
  mafsP1[i,,] <- no[i,,]/auxMP[i,]
}
mafsP2 = array(0, dim=c(dim(mafsP1)[1],dim(mafsP1)[2]*dim(mafsP1)[3]))
for (i in 1:dim(no)[1]) {
    no2[i,] = no[i,,]
    mafsP2[i,]  = mafsP1[i,,]
}
#res.pca = PCA(no2, scale.unit=TRUE, ncp=5, graph=T) 
res.pca = PCA(mafsP2, scale.unit=TRUE, ncp=5, graph=T) 
dim(mafsP2) #45 72376
# QIAGEN healthy samples ############################################################################################
pileupsQ <- list.files("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/qiagen_kit_test/201019", recursive = T, full.names = T, pattern = "bait.pileup")
pon_hg19 <- readRDS("~/genomedk/PolyA/faststorage/BACKUP/IMPROVE/call/references/200419_novaseq-xgen-sporacrc-pon.RDS") 
countsQ0 <-  piles_to_counts(files = pileupsQ, regions = pon_hg19$regions)
countsQ1 <- countsQ0[,,1:4] + countsQ0[,,6:9]
mafsQ1 = array(0, dim=c(dim(countsQ1)[1],dim(countsQ1)[2],dim(countsQ1)[3]))
auxMQ <- rowSums(countsQ1, dims = 2) 
for (i in 1:dim(countsQ1)[1]) {
  mafsQ1[i,,] <- countsQ1[i,,]/auxMQ[i,]
}
noQ = array(0, dim=c(dim(countsQ1)[1],dim(countsQ1)[2]*dim(countsQ1)[3]))
mafsQ2 = array(0, dim=c(dim(mafsQ1)[1],dim(mafsQ1)[2]*dim(mafsQ1)[3]))
for (i in 1:dim(countsQ1)[1]) {
  noQ[i,] = countsQ1[i,,]
  mafsQ2[i,]  = mafsQ1[i,,]
}
# FactoMineR: AnRPackage for Multivariate Analysis
# https://mran.microsoft.com/snapshot/2016-12-29/web/packages/FactoMineR/vignettes/FactoMineR.pdf
#res.pca = PCA(noQ, scale.unit=TRUE, ncp=5, graph=T) 
res.pca = PCA(mafsQ2, scale.unit=TRUE, ncp=5, graph=T) 
dim(mafsQ2) #24 72376
########################################################
cohort2 <- rbind(mafsP2, mafsQ2)
cohort3 <- rbind(cohort2, mafsI2)
cohort4 <- rbind(cohort3, mafsC2)
dim(cohort2) # 69 72376
res.pca = PCA(cohort2, scale.unit=TRUE, ncp=5, graph=T)
dim(cohort3) # 283 72376
res.pca = PCA(cohort3, scale.unit=TRUE, ncp=5, graph=T) 
dim(cohort4) # 373 72376
res.pca = PCA(cohort4, scale.unit=TRUE, ncp=5, graph=T) 
### IMPROVE #####
pileupsI <- list.files("~/genomedk/PolyA/faststorage/BACKUP/IMPROVE/sporacrc/N227", recursive = T, full.names = T, pattern = "bait.pileup$")
pon_hg19 <- readRDS("~/genomedk/PolyA/faststorage/BACKUP/IMPROVE/call/references/200419_novaseq-xgen-sporacrc-pon.RDS") 
countsI0 <-  piles_to_counts(files = pileupsI[1:229], regions = pon_hg19$regions)
countsI00 <- array(0, dim=c(dim(countsI0)[1]-12,dim(countsI0)[2],dim(countsI0)[3]))
countsI00[1:20,,]<-countsI0[1:20,,]
countsI00[21:30,,]<-countsI0[22:31,,] # 21
countsI00[31:45,,]<-countsI0[33:47,,] # 32
countsI00[46:49,,]<-countsI0[49:52,,] # 48
countsI00[50:56,,]<-countsI0[54:60,,] # 53
countsI00[57:113,,]<-countsI0[62:118,,] # 61
#countsI00[115,,]<-countsI0[120,,] # 119 + #120
countsI00[114:118,,]<-countsI0[123:127,,] # 121 + #122
countsI00[119:130,,]<-countsI0[129:140,,] # 128 
countsI00[131:160,,]<-countsI0[142:171,,] # 141
countsI00[161:217,,]<-countsI0[173:229,,] # 172
countsI1 <- countsI00[,,1:4] + countsI00[,,6:9]
# REMOVE #48 (N227-1981), #21 (N227-1972), #53 (N227-1985). 
mafsI1 = array(0, dim=c(dim(countsI1)[1],dim(countsI1)[2],dim(countsI1)[3]))
auxMI <- rowSums(countsI1, dims = 2) 
for (i in 1:dim(countsI1)[1]) {
  mafsI1[i,,] <- countsI1[i,,]/auxMI[i,]
}
mafsI2 = array(0, dim=c(dim(mafsI1)[1],dim(mafsI1)[2]*dim(mafsI1)[3]))
for (i in 1:dim(countsI1)[1]) {
  mafsI2[i,]  = mafsI1[i,,]
}
res.pca = PCA(mafsI2, scale.unit=TRUE, ncp=5, graph=T) 
### CRUK #####
pileupsC <- list.files("~/genomedk/PolyA/faststorage/BACKUP/CRUK/plasma/N289", recursive = T, full.names = T, pattern = "bait.pileup")
pon_hg19 <- readRDS("~/genomedk/PolyA/faststorage/BACKUP/IMPROVE/call/references/200419_novaseq-xgen-sporacrc-pon.RDS") 
countsC0 <-  piles_to_counts(files = pileupsC[1:90], regions = pon_hg19$regions)
countsC00 <- array(0, dim=c(dim(countsC0)[1]-1,dim(countsC0)[2],dim(countsC0)[3]))
countsC00[1:5,,]<-countsC0[1:5,,]
countsC00[6:89,,]<-countsC0[7:90,,]
countsC1 <- countsC00[,,1:4] + countsC00[,,6:9]
# REMOVE     #6 (N289-100)
mafsC1 = array(0, dim=c(dim(countsC1)[1],dim(countsC1)[2],dim(countsC1)[3]))
auxMC <- rowSums(countsC1, dims = 2) 
for (i in 1:dim(countsC1)[1]) {
  mafsC1[i,,] <- countsC1[i,,]/auxMC[i,]
}
mafsC2 = array(0, dim=c(dim(mafsC1)[1],dim(mafsC1)[2]*dim(mafsC1)[3]))
for (i in 1:dim(countsC1)[1]) {
  mafsC2[i,]  = mafsC1[i,,]
}
res.pca = PCA(mafsC2, scale.unit=TRUE, ncp=5, graph=T) 
###########################################################
# Hierarchical cluster analysis on a set of dissimilarities and methods for analyzing it. 
aux <- t(mafsQ2[,1:5])
hc <- hclust(aux)               
plot(hc) 
# A heatmap is another way to visualize hierarchical clustering; in the colored image, data values are transformed to color scale.
df<-scale(mafsP2[,1:5000])
col <- colorRampPalette(brewer.pal(11, "RdYlBu"))(256)
heatmap(df, scale = "none", col =  col)
# Hierarchical Clustering with P-Values via Multiscale BootstrapResampling
# https://cran.r-project.org/web/packages/pvclust/pvclust.pdf
result <- pvclust(t(mafsP2[,1:5000]), method.dist="cor", method.hclust="average", nboot=1000, parallel=TRUE)
plot(result)
