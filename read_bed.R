bed <- read.table("~/genomedk/matovanalysis/umiseq_analysis/R/NEW_METHOD_hg38_08feb2016_capture_targets.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

pos <- bed[,3] - bed[,2] + 1
sum(pos)

list <- rep(1, 18094)

bed[,1] <- bed[,3] - bed[,2] + 1

bed[1,4] = 0
bed[4,4] = 0
bed[5,4] = 0
bed[6,4] = 0
bed[15,4] = 0
bed[27,4] = 0
bed[28,4] = 0
bed[30,4] = 0
bed[31,4] = 0
bed[32,4] = 0
bed[37,4] = 0
bed[38,4] = 0
bed[42,4] = 0
bed[43,4] = 0
bed[44,4] = 0
bed[55,4] = 0
bed[56,4] = 0
bed[60,4] = 0

aux <- as.numeric(bed[,4])
aux[is.na(aux)] <- 1
bed[,4]<- aux
remove(aux)

list[1:143] = 0 # starts with a SNP
#list(1:143)<-0 # starts with a SNP

for (i in 4:60) {
  #i=3
  if (as.numeric(bed[i,4])  == 0 ) { 
    aux <- sum(bed[1:(i-1),1]) # get the current position
    a<- aux + 1
    b <- aux + bed[i,1]
    list [a:b] = 0
 print(i) 
  }
}

# 143+141+146+139+148+159+152+142+149+156+137+136+157+156+160+146+141+121
#[1] 2629
# 18094-2629
#[1] 15465