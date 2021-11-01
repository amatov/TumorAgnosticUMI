library("dplyr")

# compute W4
PON <- list.files("~/genomedk/projects/test/pon38", recursive = T, full.names = T, pattern = "txt")

PON_m = array(0, dim=c(length(PON)-1,72376/4,699))
for (i in 1:(length(PON)-1)){
  i=1
aux <- read.table(PON[i], header = F)  
auxx <- read.table(PON[i+1], header = F)  

p2 <- data.frame( aux)#PON
resC = p2 %>% summarise(errors=sum(count)) %>% filter(ref==alternate_allele)




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

aux1 <- aux[,4:702]
aux2 <- data.frame(aux1)#PON






PON_m[i,,] <- aux1
}

set.seed(8)
library(dplyr)
test <- data.frame(chr = rep("chr1", 4),
                   pos = c(1,1,2,2),
                   base = c("A", "G", "A", "G"),
                   frl1 = sample(0:10, 4),
                   frl2 = sample(0:10, 4), stringsAsFactors = F)

#Take a look at the example data
test

#Sum each columns (from column 4 to the last column) per group (group defined by chr and pos).
dplyr::group_by(test, chr, pos) %>%
  dplyr::summarise_at(.vars = names(.)[4:ncol(.)], .funs = sum)


#raw <- read.table("out.tsv", sep = "\t", header = F) #Read data
PON11[1:10, 1:10] #Take a look at a corner of the data - first 3 columns are annotation (chr, pos, base)
hist(apply(PON11[-c(1:3)], 1, sum), xlim = c(500, 3000), ylim = c(0, 10000)) #Plot, excluding 0 counts.

maxInd <- which.max(VEC)
auxVec <- PON11[maxInd, 4:702]

P1 <- PON11[list[1],4:702]
P1[116:204]
plot(seq(from = 116, to = 204, by = 1),P1[116:204])

VEC <- vector()
for (i in 1:72376){
VEC[i] <- sum(unlist(PON11[i,4:702]))
}
max(VEC) # 2149
which.max(VEC) # 67787

list <- which(VEC>0)
# 67770 67776 67777 67782 67783 67787 67790 67796 67798 67803 67806 67807 67810 67814 67819 67823 67827 67832 67836 67838 67842 67844 67847 67852 67854 67859 67862 67866 67872 67874 67877 67879 67884 67885 67891
# 67895 67899 67903 67907 67909 67911 67915 67918 67919 67924 67927 67931 67935 67938 67939 67942 67946 67950 67951 67954 67959 67961 67968 67972 67976 67977 67978 67983 67986 67991 67993 67999 68003 68007 68012
# 68016 68017 68018 68020 68023 68027 68029 68031 68035 68038 68042 68047 68049 68055 68059 68062 68063 68068 68071 68076 68080 68081 68083 68087 68090 68092 68095 68099 68104 68108 68112 68114 68116 68119 68122
# 68124 68128 68130 68132 68136 68140 68143 68146 68152 68154 68157 68164 68168 68172 68173 68180 68184 68188 68190 68195 68199 68204 68206 68209 68211 68214 68216 68220 68222 68227 68232 68234 68240 68241 68248
# 68252 68253 68258 68263 68266 68269 68274 68276 68280 68281 68288 68290 68295 68300 68302 68306 68311 68315 68317 68323 68325 68329 68336 68340 68341 68347 68351
for (i in 1:length(list)){
  VEC[list[i]]  
}

CRUKtest1 <- list.files("~/genomedk/matovanalysis/umiseq_analysis/CRUK_W3W4_CRC226", recursive = T, full.names = T, pattern = "txt")
CRUK1 <- read.table(CRUKtest1[1], header = F)  

C1 <- CRUK1[listCR[1],4:702]
C1[47:297]
plot(seq(from = 47, to = 297, by = 1),C1[47:297])

VECR <- vector()
for (i in 1:72376){
  VECR[i] <- sum(unlist(CRUK1[i,4:702]))
}

listCR <- which(VECR>0)
# 67770 67776 67777 67783 67787 67790 67792 67796 67798 67803 67807 67810 67814 67819 67822 67823 67827 67832 67834 67836 67838 67842 67844 67847 67852 67854 67858 67859 67862 67866 67872 67874 67879 67884 67885 67891
# 67895 67898 67899 67903 67907 67909 67911 67915 67919 67924 67927 67931 67935 67939 67942 67946 67951 67954 67955 67959 67961 67967 67968 67972 67976 67978 67983 67986 67988 67991 67993 67999 68003 68007 68012 68016
# 68018 68020 68023 68027 68029 68035 68038 68042 68047 68049 68055 68059 68063 68066 68068 68069 68071 68076 68077 68080 68081 68086 68087 68090 68092 68095 68099 68104 68108 68112 68114 68116 68119 68122 68128 68132
# 68136 68140 68141 68143 68146 68152 68154 68157 68164 68168 68172 68173 68180 68184 68188 68190 68195 68197 68199 68204 68205 68206 68211 68216 68220 68222 68225 68227 68232 68234 68240 68241 68248 68252 68253 68258
# 68263 68266 68269 68276 68280 68281 68288 68290 68295 68300 68302 68306 68311 68315 68317 68323 68325 68329 68336 68337 68340 68341 68347 68351
for (i in 1:length(listCR)){
  VECR[listCR[i]]  
}

max(VECR) # 2319
which.max(VECR) # 67927

length(which(CRUK1[listCR[1],4:702]>0)) #134

auxVecC <- CRUK1[which.max(VECR), 4:702]
auxVecC





