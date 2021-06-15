library("imager")
# UMICRUK
pileupsUC <- list.files("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/umiseq_paper/divergence/data/CRUK5Mb", recursive = T, full.names = T, pattern = "tsv")
umic<-vector()
nbUMIC <- length(pileupsUC)
umicD2 <- array(0,c(nbUMIC,595,499))
umic1 <-matrix(,nrow=595,ncol=499)#nrow=555,ncol=702)
j=1
for(i in 1:nbUMIC) {
  #i= 1
  print(i)
  auxFR <- read.table( pileupsUC[i], header = TRUE) # sample per sample, file per file. 
  testU <- as.integer(unlist(auxFR[,2:500]))
  umic1 <- t(matrix(testU, ncol = dim(auxFR)[1], byrow = (dim(auxFR)[2]-2)) )# convert back to matrix form
  umicD2[j,,] <- umic1
  jpeg(file="filename.jpg")
  heatmap(umicD2[j,,] )
  dev.off()

  # save as jpg
  j=j+1
  print(dim(umic1)) # 555 x 702
  if (i==1){
    umic <- umic1
  } else {
    names(umic) <- names(umic1) 
    umic <- rbind(umic , umic1)
  }
}
dim(umic) # 40460 x 499 for  68 CRUK PreOps