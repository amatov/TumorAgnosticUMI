source("~/genomedk/matovanalysis/umiseq_analysis/R/cmapply.R")
source("~/genomedk/matovanalysis/umiseq_analysis/R/image_plot.R")
#pon <- readRDS("201020_hg38-novaseq-xgen-sporacrc-pon.RDS") #hg38
pon <- readRDS("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/specs_analysis/sw_input_files/201020_hg38-novaseq-xgen-sporacrc-pon.RDS") # 

#Need the MAF in array maf
data <- pon$pon
#Turn into 18094*5*46 arrays
counts <- list()
for(pt in 1:dim(data)[[1]]) counts[[pt]] <- data[pt,,1:5] + data[pt,,6:10]
counts <- simplify2array(counts)
mafs <- list()
for(pt in 1:dim(data)[[1]]) mafs[[pt]] <- (data[pt,,1:5] + data[pt,,6:10])/rowSums(data[pt,,])
mafs <- simplify2array(mafs)

#Flag those having same as "reference"
i <- t(sapply(dimnames(mafs)[[1]], function(b) c("A", "T", "C", "G", "-") %in% b))
sum(i)
counts[i] <- NA
mafs[i] <- NA

#f1: Number of sites that have maf >= M in >= S samples
f1 <- function(a, M, S){
  sum(apply(a, c(1, 2), function(x)sum(sum(x >= M, na.rm=T) >= S)))
}

f2 <- function(a, M, S){
  #Get indeces instead of number of sites
  which(apply(a, c(1, 2), function(x)sum(sum(x >= M, na.rm=T) >= S))>0)
}

#For one comination 
S =5; M=0.01; a = mafs
indP <- f2(mafs, M, S)

#indP2 <- f2(mafs, 0.35, 1)
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

print(df)

write.table(df, "~/genomedk/matovanalysis/umiseq_analysis/blacklisting.cvs", row.names = F)

image_plot(data = log10(df[,2:ncol(df)]+1),
           plot_title = "PON sites filtering by samples and maf",
           y_title = "maf cutoff",
           breaks = seq.int(0, 4.5, 0.2),
           ylab = df[,1],
           xlab = names(df)[2:ncol(df)])
dev.off()