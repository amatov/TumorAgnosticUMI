setwd("~/genomedk/matovanalysis/umiseq_analysis/test_fragpos/")
raw <- read.table("out.txt", sep = "\t", header = F) #Read data
raw[1:10, 1:10] #Take a look at a corner of the data - first 3 columns are annotation (chr, pos, base)
hist(apply(raw[-c(1:3)], 1, sum), xlim = c(500, 3000), ylim = c(0, 10000)) #Plot, excluding 0 counts.

raw[1:20, c(1:3,150:180)]
raw[1000:1020, c(1:3,150:180)]
apply(raw[,-c(1:3)], 1, sum)


raw <- read.table("~/genomedk/projects/test/pon38/Donor294_cfdna_N140_1504_consensus.txt", sep = "\t", header = F) #Read data

raw[1:10, 4:14] #Take a look at a corner of the data - first 3 columns are annotation (chr, pos, base)
hist(apply(raw[-c(1:3)], 1, sum))#, xlim = c(500, 3000), ylim = c(0, 10000)) #Plot, excluding 0 counts.
hist(apply(raw, 1, sum))#, xlim = c(500, 3000), ylim = c(0, 10000)) #Plot, excluding 0 counts.
aux <- as.numeric(unlist(raw[,4:702]))
hist(as.numeric(unlist(raw[,4:702])))

raw[1:20, c(1:3,150:180)]
raw[1000:1020, c(1:3,150:180)]
apply(raw[,-c(1:3)], 1, sum)

pileup <- read.table('~/genomedk/PolyA/faststorage/BACKUP/IMPROVE/sporacrc/pon/201017_hg38/pon/Donor294_cfdna_N140-1504___201019153901boj/output/Donor294_cfdna_N140-1504_consensus.bait.pileup')
names(pileup)

hist(as.numeric(pileup[2:34602,5]))