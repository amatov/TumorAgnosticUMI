mahala <- function(counts0, v, list) {
  
  counts1= array(0, dim=c(dim(counts0)[1],sum(list),dim(counts0)[3]))
  for (i in 1:dim(counts)[1]) {
    p2 <- data.frame(counts0[i,,])#CRUK
    p1 <- p2 [list == 1, ] 
    counts1[i,,] <- data.matrix(p1)
  }
  counts <- counts1[,,1:4] + counts1[,,6:9]
  
  mafs = array(0, dim=c(dim(counts)[1],dim(counts)[2],dim(counts)[3]))
  auxM <- rowSums(counts, dims = 2) 
  for (i in 1:dim(counts)[1]) {
    mafs[i,,] <- counts[i,,]  /auxM[i,]
  }
  w1 = array(0, dim=c(dim(counts)[1],dim(counts)[2],dim(counts)[3]))
  sc <- vector()
  for (i in 1:dim(mafs)[1]) {
    print(i)
    #w1[i,,] <- mafs[i,,]#/v/sum(list)  # the patient MAFs devided by the variance of the position based on PON'
    sc[i]<-sum(w1[i,,])
  }
  return(sc) # mafs, counts
}