pileupsQ <- list.files("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/qiagen_kit_test/201019", recursive = T, full.names = T, pattern = "bait.pileup")
countsQ00 <-  piles_to_counts(files = pileupsQ, 
                             regions = pon_hg19$regions)
countsQ = array(0, dim=c(dim(countsQ00)[1]-2,dim(countsQ00)[2],dim(countsQ00)[3]))
countsQ[1:16,,] <- countsQ00[1:16,,]
countsQ[17:22,,]<-countsQ00[19:24,,]
counts <- countsQ
############################
pon_obj2 <- readRDS("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/umiseq_paper/reference/201217_hg38-novaseq-xgen-sporacrc-pon.RDS") # 46
#pon_obj2 <- readRDS("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/specs_analysis/sw_input_files/201020_hg38-novaseq-xgen-sporacrc-pon.RDS") # 45
pon_counts <- pon_obj2[["pon"]]
no = array(0, dim=c(dim(pon_counts)[1]-1,dim(pon_counts)[2],dim(pon_counts)[3]))
no[1:27,,] <- pon_counts[1:27,,]
no[28:44,,]<-pon_counts[29:45,,]
counts <- no
##############################
#Turn into 18094*5*46 arrays
library(dplyr)
mafs <- abind::abind(
  lapply(1:dim(counts)[[1]], function(s)(counts[s,,1:5]+counts[s,,6:10])/rowSums(counts[s,,])),
  along = 0
)
#Indeces of the IDSNPs
coordinates <- pon$coordinates
regions <- read.table("~/genomedk/PolyA/faststorage/BACKUP/IMPROVE/sporacrc/sporacrcv1_bed/NEW_METHOD_hg19_08feb2016_capture_targets.bed")
#This is slow ~2 min O ~ n² (but works)
i <-
  sapply(1:nrow(coordinates), function(x){
    any(sapply(1:nrow(regions), function(y) {
      grepl("IDSNP", regions[y, 4]) &
        coordinates[x, 1] == regions[y, 1] &
        coordinates[x, 2] >= regions[y, 2] &
        coordinates[x, 2] <= regions[y, 3]
    }
    )
    )
  }
  )

sum(i) == sum((regions$V3 - regions$V2 + 1)[grepl("IDSNP", regions$V4)]) #TRUE

#Indeces of reference
j <- t(sapply(coordinates$ref, function(b) c("A", "T", "C", "G", "-") %in% b))

#All possible SNV sitemuts on panel (including DELs) as types
sitemuts_panel <- t(apply(pon$coordinates, 1, function(x){
  paste(x[3], c("A", "T", "C", "G", "Del" ), sep = ">")}))
#pon$genome
subs <- as.vector(unique(sitemuts_panel))
subs <- subs[!sub("(.+)>(.+)", "\\1", subs) == (sub("(.+)>(.+)", "\\2", subs))]

data <- list()
for (s in 1:length(subs)) {
  # where subs <- as.vector(unique(sitemuts_panel))
  i0 <- which(sitemuts_panel == subs[s])
  i0 <- setdiff(i0, union(which(j), which(i)))
  data[[subs[s]]] <-  as.vector(apply(mafs, 1, function(x) x[ i0 ]))
}

pdata <- sapply(data, `length<-`, max(lengths(data)))
pdata <- as.data.frame(pdata)

pdata_m <- 
  reshape2::melt(pdata, measure.var = colnames(pdata)) %>% 
  dplyr::filter(value < 0.1, 
                value > 0)

boxplot(value ~ variable, data = pdata_m, log = "y", ylab = "MAF", xlab = "", 
        main = "MAFs of core positions in all PONs ]0; 0.1[")
lbl <- group_by(pdata_m, variable) %>% summarise(n = n())
mtext(paste0(trimws(format(round(lbl[,2, drop = T]/1000, 1), nsmall = 1)), "k" ),
      side = 1,
      at = 1:nrow(lbl), 
      line = -1)

pdata_m <- 
  reshape2::melt(pdata, measure.var = colnames(pdata)) %>% 
  dplyr::filter(value < 0.1)
# CHANGE TO VIOLIN PLOT 
boxplot(value ~ variable, data = pdata_m, ylab = "VAF", xlab = "", 
        main = "MAFs of core positions in all PONs [0; 0.1[")



#Same but use mean across PON
data <- list()
for (s in 1:length(subs)) {
  i0 <- which(sitemuts_panel == subs[s])
  i0 <- setdiff(i0, union(which(j), which(i)))
  data[[subs[s]]] <-  rowMeans(apply(mafs, 1, function(x) x[ i0 ]))
}


#list of vectors to matrix filling out with NA
pdata <- sapply(data, `length<-`, max(lengths(data)))

pdata <- as.data.frame(pdata)
pdata_m <- 
  reshape2::melt(pdata, measure.var = colnames(pdata)) %>% 
  dplyr::filter(value < 0.1, 
                value > 0)
lbl <- group_by(pdata_m, variable) %>% summarise(n = n())

boxplot(value ~ variable, data = pdata_m, log = "y", ylab = "VAF", xlab = "", 
        main = "Mean VAFs of core positions across 22 QIAGEN samples ]0; 0.1[")
mtext(paste0(trimws(format(round(lbl[,2, drop = T]/1000, 1), nsmall = 1)), "k" ),
      side = 1,
      at = 1:nrow(lbl), 
      line = -1)

pdata_m <- 
  reshape2::melt(pdata, measure.var = colnames(pdata)) %>% 
  dplyr::filter(value < 0.1)
# CHANGE TO VIOLIN PLOT 
boxplot(value ~ variable, data = pdata_m, ylab = "VAF", xlab = "", 
        main = "Mean VAFs of core positions across 22 QIAGEN samples [0; 0.1[")

pdata1 <- as.data.frame(pdata)
p <- ggplot(pdata, aes(x=colnames(pdata), y=rownames(pdata))) + geom_violin()