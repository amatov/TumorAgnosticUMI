library("dplyr")
# PON mutations and variability ###################################################################################
pon_obj2 <- readRDS("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/specs_analysis/sw_input_files/201020_hg38-novaseq-xgen-sporacrc-pon.RDS") # 
sitemut <- t(apply(pon_obj2$coordinates, 1, function(x){
  paste( paste0(trimws(x[1]), ":", trimws(x[2]), "_"),
         paste(x[3] ,c("A", "T", "C", "G"), sep = "/"),
         sep = "")}))
cou <- pon_obj2[["pon"]]
co <- cou[,,1:4]+cou[,,6:9]
vo<-apply(co,2:3,var) # variance of the counts of each position based on PON
v0 <- min(vo[vo>0])/10000000 # for counts w zero variance, we replace w a very small value
vo1<- vo
vo1[vo==0]=v0
# COSMIC prior ##########################################################################################
prior1 <- readRDS("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/specs_analysis/sw_input_files/180903_prior.RDS")
prior11 <- prior1[,1:4]
# IMPROVE plasma data ##################################################################################################
pileupsIw <- list.files("~/genomedk/PolyA/faststorage/BACKUP/IMPROVE/sporacrc/N227", recursive = T, full.names = T, pattern = "bait.pileup")
source("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/specs_analysis/sw_input_files/duplex_tools.R")
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