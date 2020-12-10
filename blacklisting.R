source("~/genomedk/matovanalysis/umiseq_analysis/R/read_bed.R") #1/0 list
source("~/genomedk/matovanalysis/umiseq_analysis/R/cmapply.R")
source("~/genomedk/matovanalysis/umiseq_analysis/R/image_plot.R")
#pon <- readRDS("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/specs_analysis/sw_input_files/201020_hg38-novaseq-xgen-sporacrc-pon.RDS") # 
pon_obj2 <- readRDS("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/specs_analysis/sw_input_files/201020_hg38-novaseq-xgen-sporacrc-pon.RDS") # 
pon_counts <- pon_obj2[["pon"]]
inRef <- t(sapply(dimnames(pon_counts)[[2]], function(b) c("A", "T", "C", "G") %in% b))
sum(inRef)
#pon_counts[inRef] <- NA
#no1 = array(0, dim=c(dim(pon_counts)[1],sum(list),dim(pon_counts)[3]))
no <- pon_counts[,,1:4]+pon_counts[,,6:9]
mafsP1 = array(0, dim=c(dim(no)[1],dim(no)[2],dim(no)[3]))
auxMP <- rowSums(no, dims = 2) 
for (i in 1:dim(no)[1]) {
  mafsP1[i,,] <- no[i,,]  /auxMP[i,]
}
#mafsP2 = array(0, dim=c(dim(no)[1],sum(list),dim(no)[3]))# core panel
mafsP2 = array(0, dim=c(dim(no)[1],dim(no)[2],dim(no)[3]))# full list
for (i in 1:dim(no)[1]) {
  #i=1
  auxRP <- mafsP1[i,,]
  auxRP[inRef] <- NA
  p2 <- data.frame(auxRP)#PON
  p1 <- p2 [list == 1, ] 
  mafsP2[i,,] <- data.matrix(p1)
}

f2 <- function(a, M, S){
  #Get indexes instead of number of sites
  which(apply(a, c(1, 2), function(x)sum(sum(x >= M, na.rm=T) >= S))>0)
  #which(apply(a, c(2, 3), function(x)sum(sum(x >= M, na.rm=T) >= S))>0)
}
f3 <- function(a, M, S){
  #Get indexes instead of number of sites
  which(apply(a, c(2, 3), function(x)sum(sum(x >= M, na.rm=T) >= S))>0)
  #which(apply(a, c(2, 3), function(x)sum(sum(x >= M, na.rm=T) >= S))>0)
}
#apply(mafs[1:1000,,], c(1,2), function(x)str(x))
#apply(mafsP2[,1:1000,], c(2,3), function(x)str(x))
#For one combination 
S =5; M=0.01;  
f3(mafsP2, M, S)
f2(mafs, M, S)# ref included.

indP <- f2(mafs, 0.01, 10)
#setdiff(indP, indP2)
#which(  , arr.ind = T )


mafs[ which( mafs[,,1] > 0.1, arr.ind = T)[,1], ,1]

mafs[,,1]

sort(mafs[,,1][indP], decreasing = TRUE)

#tmp <- apply(a, c(1, 2), function(x)sum(sum(x >= M, na.rm=T) >= S))
#which(tmp > 0)

#Use cmapply do do all combinations of arguments M and S
r <- cmapply(f1,
             M = c(0.1, 0.5, 1, 5,10, 35)/100,#, 65, 85, 100)/100,
             S = c(3,4,5,6),#1, 10, 20, 30, 40, 43, 46),
             MoreArgs = list(a = mafs))

df <- reshape2::dcast(r, M ~ S)

#MADS
# mean(mafs[mafs<0.1 & mafs >0], na.rm = T)
mafsO <- mafs
mean(mafsO[ mafsO<=0.35 ], na.rm = T) # 3.044935e-05
for (i in 1:45) {
  mafs[,,i][indP]<-NA
}
mean(mafs[  mafs<=0.35 ], na.rm = T) # 2.322167e-05
er <- vector()
erO <- vector()
for (i in 1:45) {
  auxB <- mafs[,,i] 
  er[i] <- mean(auxB[auxB<=0.35], na.rm=T)
  auxBO <- mafsO[,,i] 
  erO[i] <- mean(auxBO[auxBO<=0.35], na.rm=T)
}

length(mafs) 


print(df)

write.table(df, "~/genomedk/matovanalysis/umiseq_analysis/blacklisting.cvs", row.names = F)

image_plot(data = log10(df[,2:ncol(df)]+1),
           plot_title = "PON sites filtering by samples and maf",
           y_title = "maf cutoff",
           breaks = seq.int(0, 4.5, 0.2),
           ylab = df[,1],
           xlab = names(df)[2:ncol(df)])
dev.off()

pon <- readRDS("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/specs_analysis/sw_input_files/201020_hg38-novaseq-xgen-sporacrc-pon.RDS") # 
data <- pon$pon
#Turn into 18094*5*46 arrays
counts <- list()
for(pt in 1:dim(data)[[1]]) counts[[pt]] <- data[pt,,1:5] + data[pt,,6:10]
counts <- simplify2array(counts)
mafs <- list()
for(pt in 1:dim(data)[[1]]) mafs[[pt]] <- (data[pt,,1:5] + data[pt,,6:10])/rowSums(data[pt,,])
mafs <- simplify2array(mafs)

#Flag those having same as "reference"
inM <- t(sapply(dimnames(mafs)[[1]], function(b) c("A", "T", "C", "G", "-") %in% b))
sum(inM)
counts[inM] <- NA
mafs[inM] <- NA