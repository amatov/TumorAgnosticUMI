countsW <- readRDS("~/genomedk/matovanalysis/umiseq_analysis/R/cruki.RDS") # # REPLACE WITH YOUR COPY

# 45 Subjects of the Control Panel of Normal PON ####################################################################
pon_obj2 <- readRDS("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/umiseq_paper/reference/201217_hg38-novaseq-xgen-sporacrc-pon.RDS") # 46
pon_counts <- pon_obj2[["pon"]]
no0 = array(0, dim=c(dim(pon_counts)[1]-1,dim(pon_counts)[2],dim(pon_counts)[3]))
no0[1:27,,] <- pon_counts[1:27,,]
no0[28:45,,]<-pon_counts[29:46,,]
no <- no0[,,1:4]+no0[,,6:9]
mafsP1 = array(0, dim=c(dim(no)[1],dim(no)[2],dim(no)[3]))
auxMP <- rowSums(no, dims = 2) 
for (i in 1:dim(no)[1]) {
  mafsP1[i,,] <- no[i,,]  /auxMP[i,]
}
############################# Variance of the PON at each of the 62k positions #######################################3
v<-apply(mafsP1,2:3,sd) # variance of the counts of each position based on PON
v0 <- min(v[v>0])/100000#00 # for counts w zero variance, we replace w a very small value
v1<- v
v1[v==0]=v0
# PON mutations and variability ###################################################################################
sitemut <- t(apply(pon_obj2$coordinates, 1, function(x){
  paste( paste0(trimws(x[1]), ":", trimws(x[2]), "_"),
         paste(x[3] ,c("A", "T", "C", "G"), sep = "/"),
         sep = "")}))

# cruk plasma data ##################################################################################################
countsC00 <- readRDS("~/genomedk/matovanalysis/umiseq_analysis/R/cruk-counts.RDS") # REPLACE WITH YOUR COPY
countsC001 <- countsC00[,,1:4] + countsC00[,,6:9]
mafsC = array(0, dim=c(dim(countsC001)[1],dim(countsC001)[2],dim(countsC001)[3]))
auxC <- rowSums(countsC001, dims = 2) 
for (i in 1:dim(countsC001)[1]) {
  mafsC[i,,] <- countsC001[i,,]  /auxC[i,] 
}
#this is from your code to generate cancer_SNPs:
#info <- read.table("~/genomedk/matovanalysis/umiseq_analysis/R/cruk-plasma-info.lst", header = T, stringsAsFactors = F) # REPLACE W YOUR COPY
#tmp <- info[ nchar(sub(p, "\\1\\2", info$sitemut_hg38)) == 0, ]
#cancer_SNPs <- tmp[!is.na(tmp$sitemut_hg38), ]

cancer_SNPs$i # 45 unique
#  6  6  6  7  8  9  9 10 10 10 10 11 11 11 12 12 12 12 18 18 18 18 18 21 21 21 22 22 22 22 23 23 24 24 24 24 24 25 25 26 26 26 26 27 27
# 27 28 28 31 31 31 32 32 32 32 33 33 33 35 35 37 37 37 37 38 38 39 39 41 43 43 43 44 44 44 45 46 46 48 48 49 49 50 50 50 51 51 54 54 54
# 54 55 55 56 56 58 58 58 60 60 60 62 63 64 65 66 66 67 67 67 68 68 68 69
cancer_list <- unique(cancer_SNPs$i)
# 6  7  8  9 10 11 12 18 21 22 23 24 25 26 27 28 31 32 33 35 37 38 39 41 43 44 45 46 48 49 50 51 54 55 56 58 60 62 63 64 65 66 67 68 69

maha <- vector()
j=1
for (i in cancer_list){
  maha[j] <- sum(countsW[i,,]* mafsC[i,,]/v1,  na.rm=T) / sum(countsW[i,,])
} 
print(maha)
 