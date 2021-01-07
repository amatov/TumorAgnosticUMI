# Method 1: bbb_ calls and describe founder samples by unfiltered bbb_ on pileups -------------
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()
#install.packages()
library(matrixStats)
library(reshape2)
#setwd ('G:\\PolyA/faststorage/BACKUP/N140_Targeting/specs/specs_analysis')
setwd ('~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/specs_analysis')
library(ROCR)
source("U:\\Documents/R/utility_functions-master/recoder.R")
source("U:\\Documents/R/utility_functions-master/auc.R")
source("U:\\Documents/R/utility_functions-master/scaler.R")
source("U:\\Documents/R/utility_functions-master/confusion_plot.R")
source("G:\\PolyA/faststorage/BACKUP/N140_Targeting/specs/specs_analysis/sw_piles.R")
source("sw_input_files/duplex_tools.R")
library("dplyr")
library("tidyr")
library("ggplot2")
library(qlcMatrix) 
library(FactoMineR)
source("~/genomedk/matovanalysis/umiseq_analysis/R/read_bed.R")
############################################################################################################################
pon_hg19 <- readRDS("~/genomedk/PolyA/faststorage/BACKUP/IMPROVE/call/references/200419_novaseq-xgen-sporacrc-pon.RDS") # 
str(pon_hg19)
# 45 Subjects of the Control Panel of Normal PON
pon_obj2 <- readRDS("sw_input_files/201020_hg38-novaseq-xgen-sporacrc-pon.RDS") # 
pon_counts <- pon_obj2[["pon"]]
str(pon_obj2)
no = pon_counts[,,1:4]+pon_counts[,,6:9] # no = pon_counts[,,1:5]+pon_counts[,,6:10]
# QIAGEN control 24 samples HG19 ################################################################################################
pileupsQ <- list.files("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/qiagen_kit_test/201019", recursive = T, full.names = T, pattern = "bait.pileup")
all(file.exists(pileupsQ))
countsQ <-  piles_to_counts(files = pileupsQ, 
                           regions = pon_hg19$regions)
countsQ1 <- countsQ[,,1:4] + countsQ[,,6:9]# counts for the 24 Qiagen samples
mQ1 <- cQ1/rowSums(cQ1)
mQ11 <- mQ1[list == 1]
mQ111 <- mQ11 [mQ11 <= 0.1]
#mQ1111 <- mQ111 [mQ111 > 0 ]
plot(sort(mQ111 [1:100], decreasing = TRUE))

plot(sort(mQ11 [mQ11 <1]),decreasing = TRUE)
plot(sort(m1), decreasing = TRUE)

cQ1 <- countsQ1[1,,] 
cQ11 <- matrix(cQ1[list == 1],nrow=sum(list),ncol=4)
plot(sort(rowSums(cQ11), decreasing = TRUE)) # Coverage per position for 1st sample

covQ <- rowSums(countsQ1, dims = 2)
meaQ<-apply(covQ,2,mean) # variance of each position based on PON
plot(sort(meaQ, decreasing = TRUE))
covQ1 <- t(rbind(covQ,meaQ))
covQ2 <- covQ1[order(covQ1[,25],decreasing=TRUE),]
plot(covQ2[,25])# the sorted mean coverage for each position
covQ3 <- covQ2[,1:24]# the sorted coverage for each position in the 24 samples
#########Gray band (variance) ############################
means = apply(covQ3,1,mean)
sds = apply(covQ3,1,sd)
df = data.frame(idx = 1:nrow(covQ3), mean = means, sd = sds)
#df2 = data.frame(sort(rowSums(cQ11), decreasing = TRUE))  
df2 = data.frame(sort(means, decreasing = TRUE))  
ggplot(df, aes(x = idx, y = mean)) +geom_line() + geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd) , fill = "grey40",alpha=.5) 
ggplot(df, aes(x = idx, y = mean)) + geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd) ) 
ggplot(df, aes(x = idx, y = mean)) + geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd) , fill = "grey70") + geom_step(data = df2)
# SHOW CI ###############################################
means <- covQ2[,25]
stdev <- sqrt(apply(covQ3,1,var))
n     <- dim(covQ3)[1]
ciw   <- qt(0.975, n) * stdev / sqrt(n)
plotCI(x=means, uiw=ciw, col="black", barcol="blue", lwd=1)
########################################
mafsQ1 <- countsQ1/rowSums(countsQ1)
covQ1 <- rowSums(countsQ1)
mf <- mafsQ1[1,,]
mfl <- matrix(mf[list == 1],nrow=sum(list),ncol=4)
plot(sort(mfl[mfl>0], decreasing = TRUE))
plot(sort(mf[mf>0], decreasing = TRUE))

# CRUK control 8 samples #################################################################################################
pileupsC <- list.files("~/genomedk/PolyA/faststorage/BACKUP/CRUK/plasma/N289", recursive = T, full.names = T, pattern = "bait.pileup")
countsC <-  piles_to_counts(files = pileupsC[1:8], 
                            regions = pon_hg19$regions)
countsC1 <- countsC[,,1:4] + countsC[,,6:9]# counts for the 24 Qiagen samples

# Eight Panels of Dilution Series###########################################################################################
dat0 <- readRDS("sw_output_files/2020-10-23-145546_sw-output.RDS") # DS SW output hg38
pileupsD <- unlist(attributes(dat0))
pileupsD <- sub("/faststorage/project/PolyA/BACKUP", "~/genomedk/PolyA/faststorage/BACKUP", pileupsD)
all(file.exists(pileupsD))
counts <-  piles_to_counts(files = pileupsD, 
                          regions = pon_obj2$regions) 
counts1 <- counts[,,1:4] + counts[,,6:9]# counts for the 8 DS panels

# noise calculations
r1 <-rownames(no[1,,])
p2 <- data.frame(ref=r1, no[1,,])#PON
#p2 <- data.frame(ref=r1, counts1[1,,]) #DS 
#countsQ1[1,,][indP]<-NA # blacklist
#p2 <- data.frame(ref=r1, countsQ1[1,,])#QIAGEN
#p2 <- data.frame(ref=r1, countsC1[1,,])#CRUK
p1 <- p2 [list == 1, ] 
 
res1 = p1 %>%
  mutate(Allele_reads = A + T + C + G
  ) %>%
  filter(Allele_reads>0) %>%
  mutate(A_freq = A/Allele_reads,
         T_freq = T/Allele_reads,
         C_freq = C/Allele_reads,
         G_freq = G/Allele_reads) %>%
  data.frame(.)

resP = p1 %>%
  pivot_longer(cols = matches("^[ATCG]"),
               names_to = c("alternate_allele"),
               values_to = "count", 
               names_pattern = "([ATCG])")

resCH = resP %>% filter(ref!=alternate_allele)
#ggplot(resP %>% filter(ref!=alternate_allele), aes(x = count)) + geom_histogram()

#resP$sampleID = 1
resC = resP  %>% group_by(ref, alternate_allele) %>% summarise(errors=sum(count)) %>% filter(ref==alternate_allele)
C <- sum(resC$errors)
resE = resP  %>% group_by(ref, alternate_allele) %>% summarise(errors=sum(count)) %>% filter(ref!=alternate_allele)
E <- sum(resE$errors)
Er <- E/(C+E)*100
#resE = resP  %>% filter(ref!=alternate_allele) %>% group_by(ref, alternate_allele) %>% summarise(errors=sum(count)) 

res = resP  %>% group_by(ref, alternate_allele) %>% summarise(errors=sum(count)) %>% filter(ref!=alternate_allele) %>% ungroup() %>% mutate(error_freq = errors /sum(errors))
res$sampleID = 1
res$nonrefcounts = E

#samples= c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,19,20,21,22,23,24)
#samples= c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24, 25, 26, 27, 29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45)

for (i in 2:45) {
  print(i)
  p2 <- data.frame(ref=r1, no[i,,])#PON
  #p2 <- data.frame(ref=r1, counts1[i,,])#DS
  #countsQ1[i,,][indP]<-NA # blacklist
  #p2 <- data.frame(ref=r1, countsQ1[i,,])#QIAGEN
  #p2 <- data.frame(ref=r1, countsC1[i,,])#CRUK
  p1 <- p2 [list == 1, ] 
  
res1 = p1 %>%
  mutate(Allele_reads = A + T + C + G
  ) %>%
  filter(Allele_reads>0) %>%
  mutate(A_freq = A/Allele_reads,
         T_freq = T/Allele_reads,
         C_freq = C/Allele_reads,
         G_freq = G/Allele_reads) %>%
  data.frame(.)

resP = p1 %>%
  pivot_longer(cols = matches("^[ATCG]"),
               names_to = c("alternate_allele"),
               values_to = "count", 
               names_pattern = "([ATCG])")

resCH1 = resP %>% filter(ref!=alternate_allele)
resCH <- bind_rows(resCH , resCH1 ) 
#resP$sampleID = i
resC1 = resP  %>% group_by(ref, alternate_allele) %>% summarise(errors=sum(count)) %>% filter(ref==alternate_allele)
C1 <- sum(resC1$errors)
resE1 = resP  %>% group_by(ref, alternate_allele) %>% summarise(errors=sum(count)) %>% filter(ref!=alternate_allele)
E1 <- sum(resE1$errors)
Er1 <- E1/(C1+E1)*100
Er <- append(Er, Er1)

res1 = resP  %>% group_by(ref, alternate_allele) %>% summarise(errors=sum(count)) %>% filter(ref!=alternate_allele) %>% ungroup() %>% mutate(error_freq = errors /sum(errors))
res1$sampleID = i
res1$nonrefcounts = E1
#res <- dplyr::full_join(res, res1, by = c("ref", "alternate_allele")) 
  
#res <- bind_rows(res, res1, by = c("ref", "alternate_allele"))  
res <- bind_rows(res, res1) 
resE <- bind_rows(resE, resE1) 
}
plot(Er,pch = 19, ylim=c(0.02,0.08))
ErI <- Er
#ggplot(res, aes(x = sampleID, y = error_freq)) + geom_bar(position="fill", stat="identity")

Res_plot = res %>% mutate(mutation = paste0(ref,">",alternate_allele))
#df_total = Res_plot  %>% group_by(sampleID) %>% summarize(total_count = sum(count))
ggplot(Res_plot, aes(x = as.factor(sampleID), y = error_freq, fill=mutation)) + geom_bar(position="fill", stat="identity", color = "black") #+ geom_text(aes(x = sampleID, y = 0, label=nonrefcounts, angle = 90,size = 2))

# absolute errors
ggplot(Res_plot, aes(x = as.factor(sampleID), y = errors, fill=mutation)) + geom_bar( stat="identity", color = "black")

# check ref. for VAF>0.5
resCHr = resP %>% filter(ref==alternate_allele)
test <- resP %>% mutate(new_feature = nonrefcounts > ref)
refCount <- resCHr$count # ref vector 15,465
nrfCount <- resCH1$count
nrmCount <-matrix(nrfCount,nrow=15465,ncol=3)
nrM <- rowMax(nrmCount) # non ref vector 15,465
difV <- refCount - nrM
which(difV<0) # PON45:  3085  4688  7245 12725

#ggplot(Res_plot, aes(x = sampleID, y = error_freq, fill=mutation)) + geom_bar(position="fill", stat="identity", color = "black") + theme(axis.text.x = element_text(angle = 90)) 
#+ labs(x="", y="(%)")

#ggplot(Res_plot, aes(x = sampleID, y = error_freq, fill=mutation)) + geom_bar(position="stack", stat="identity", color = "black")
#ggplot(res, aes(x = sampleID, y = error_freq, col =mutation)) +  geom_bar(position="fill", stat="identity", color = "black")

# Look for CH based on unimodal/bimodal transition histogram
resCH2 <- resCH %>% filter(count>0)
ggplot(resCH2, aes(x = count)) + geom_histogram()
CH <- unlist(resCH2[,3])
#ggplot(Res_plot, aes(x = count)) + geom_histogram() + facet_wrap(vars(mutations_type))
#ggplot(data, aes(x = count, group = mutation_type, color = mutation_type)) + density()

# nice colors but no separation contours
ggplot(res, aes(x = sampleID, y = error_freq, fill = interaction(ref, alternate_allele))) + geom_bar(position="fill", stat="identity")

# black with contours around each transition
ggplot(res, aes(x = sampleID, y = error_freq, col = interaction(ref, alternate_allele))) + geom_bar(position="fill", stat="identity")

#ggplot(res, aes(x = sampleID, y = error_freq, color = interaction(ref, alternate_allele)) + geom_bar(position="fill"))


setwd ('U:\\Documents/R')
library(xlsx)
#write.xlsx(x, file, sheetName="Sheet1")
#write.csv(Your DataFrame,"Path where you'd like to export the DataFrame\\File Name.csv", row.names = FALSE)
write.csv(res, "U:\\Documents/R/PON_transitions.csv",na="NA",row.names=TRUE)


 