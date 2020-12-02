library(FactoMineR)
# 45 Subjects of the Control Panel of Normal PON ####################################################################
pon_obj2 <- readRDS("sw_input_files/201020_hg38-novaseq-xgen-sporacrc-pon.RDS") # 
no <- pon_obj2[["pon"]] # counts
no2 = array(0, dim=c(dim(no)[1],dim(no)[2]*dim(no)[3]))
mafsP2 = array(0, dim=c(dim(mafsP1)[1],dim(mafsP1)[2]*dim(mafsP1)[3]))
for (i in 1:45) {
    no2[i,] = no[i,,]
    mafsP2[i,]  = mafsP1[i,,]
}
res.pca = PCA(no2, scale.unit=TRUE, ncp=5, graph=T) 
res.pca = PCA(mafsP2, scale.unit=TRUE, ncp=5, graph=T) 
# QIAGEN healthy samples ############################################################################################
pileupsQ <- list.files("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/qiagen_kit_test/201019", recursive = T, full.names = T, pattern = "bait.pileup")
countsQ0 <-  piles_to_counts(files = pileupsQ, 
                             regions = pon_hg19$regions)
countsQ= array(0, dim=c(dim(countsQ0)[1],sum(list),dim(countsQ0)[3]))
noQ = array(0, dim=c(dim(countsQ0)[1],dim(countsQ0)[2]*dim(countsQ0)[3]))
mafsQ2 = array(0, dim=c(dim(mafsQ1)[1],dim(mafsQ1)[2]*dim(mafsQ1)[3]))
for (i in 1:24) {
  noQ[i,] = countsQ0[i,,]
  mafsQ2[i,]  = mafsQ1[i,,]
}
res.pca = PCA(noQ, scale.unit=TRUE, ncp=5, graph=T) 
res.pca = PCA(mafsQ2, scale.unit=TRUE, ncp=5, graph=T) 
########################################################
 

