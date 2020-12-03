source("sw_input_files/duplex_tools.R")
pon_obj2 <- readRDS("sw_input_files/201020_hg38-novaseq-xgen-sporacrc-pon.RDS") # 
cou <- pon_obj2[["pon"]]
co <- cou[,,1:4]+cou[,,6:9]
vo<-apply(co,2:3,var) # variance of the counts of each position based on PON
v0 <- min(vo[vo>0])/10000000 # for counts w zero variance, we replace w a very small value
vo1<- vo
vo1[vo==0]=v0
prior1 <- readRDS("sw_input_files/180903_prior.RDS")
prior11 <- prior1[,1:4]
###################
maW = array(0, dim=c(dim(co)[1],dim(co)[2],dim(co)[3]))
rS <- rowSums(co, dims = 2) 
for (i in 1:dim(co)[1]) {
  #i = 1
  maW[i,,] <- co[i,,]/rS[i,]
}
#All possible SNV sitemuts on panel##################################################################
sitemut <- t(apply(pon_obj2$coordinates, 1, function(x){
  paste( paste0(trimws(x[1]), ":", trimws(x[2]), "_"),
         paste(x[3] ,c("A", "T", "C", "G"), sep = "/"),
         sep = "")}))
# IMPROVE WES data ###################################################################################
wes <- read.table("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/umiseq_paper/data/201123_wes-spora-mutations-improve.csv", header = TRUE) # HG38
#wes <- read.table("~/genomedk/matovanalysis/umiseq_analysis/201123_wes-spora-mutations-improve.csv", header = TRUE) # HG38
pts <- unique(wes$pt_id)
#sapply(as.character(pts), function(x) grep(x, pileupsIw[preop_i])) # 141 of 179
#wes_id <- intersect(ItW$pt_id, pts)
ItW <- read.table("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/umiseq_paper/data/IMPROVEptList",header = TRUE)
#ItW <- read.table("~/genomedk/matovanalysis/umiseq_analysis/IMPROVEptList",header = TRUE)
wes_id <- ItW$pt_id %in% pts
#ItW$index <- sapply(as.character(ItW$library_id), function(x) grep(x, pileupsIw)) # 141 of 179
preI <- (ItW$op_time_cat == -1)&wes_id
preop_i <- ItW[preI, "index"] # 57
posop14_i <- ItW[ (ItW$op_time_cat == 2)&wes_id, "index"] # 42
posop30_i <- ItW[ (ItW$op_time_cat == 30)&wes_id, "index"] # 42

#mutations for the individual pt
wesP <- wes[wes==pts[1],]
m1 <- wesP$sitemut_hg38[1]
#chr17:7675235_T/C
#chr3:179218294_G/A
#chr7:140753336_A/T
sum(sitemut == "chr17:7675235_T/C") # 1
sum(sitemut == m1) # 1
which(sitemut == m1)# 50720
sitemut[which(sitemut == m1)] # "chr17:7675235_T/C"
###################################################################################################
pileupsIw <- list.files("~/genomedk/PolyA/faststorage/BACKUP/IMPROVE/sporacrc/N227", recursive = T, full.names = T, pattern = "bait.pileup")
countsIw <-  piles_to_counts(files = pileupsIw[preop_i], regions = pon_hg19$regions) # PREOP, POSTOP14, POSTOP30
countsW <- countsIw[,,1:4] + countsIw[,,6:9]
mafsW= array(0, dim=c(dim(countsW)[1],dim(countsW)[2],dim(countsW)[3]))
auxM <- rowSums(countsW, dims = 2) 
for (i in 1:dim(countsW)[1]) {
  mafsW[i,,] <- countsW[i,,]  /auxM[i,]
}
ww = array(0, dim=c(dim(co)[1],dim(co)[2],dim(co)[3]))
sw<- vector()
#m<-wesP$sitemut_hg38
wQ = array(0, dim=c(dim(countsW)[1],dim(countsW)[2],dim(countsW)[3]))
scQ <- vector()
counter <- 1
for (i in 1:dim(mafsW)[1]) {
  i=1
  j= 95
  for (j in 1:length(pts)) {  
    if (grepl(as.character(pts[j]), as.character(pileupsIw[preop_i][i]))) {# TRUE
      #print(i)
      #print(j)
      print(counter)
      counter<- counter+1
      wesP <- wes[wes==pts[j],]
      mu <- wesP$sitemut_hg38
      for (k in 1:length(mu)) {
        auxW <- mafsW[i,,]
        iW <- sitemut == mu[k] # sum(iW) is 1
        wQ[k] <- auxW[iW]/vo1[iW]*prior11[iW]
      }
      scQ[j] <- sum(wQ)/length(mu)
    }
}
}
for (i in 1:length(m)) {
  i = 3
  # index of mutations in the list
  which(sitemut == m[i])
  ww[i,,] <- mafsQ1[i,,][mafsQ1[i,,]<= 0.1]/v1/sum(mafsQ1[i,,]<= 0.1)*prior2
  sw[i]<-sum(ww[i,,])
}
#######################################################################
r1<- -3
r2 <- 14
plot(log2(sc1Q), ylim=range(c(r1,r2)), col="green", pch = 17)
par(new = TRUE)
plot(log2(sc11Q), ylim=range(c(r1,r2)), col="red", pch = 19)
par(new = TRUE)
plot(log2(sc111Q), ylim=range(c(r1,r2)), col="blue", pch = 18)