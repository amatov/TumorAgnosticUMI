### Intersect WES mutations with plasma and calculate score ###-----------------
library("dplyr")

wes <- read.table("~/genomedk/N140_Targeting/specs/umiseq_paper/data/201123_wes-spora-mutations-improve.csv", header = TRUE) # HG38
plasma <- read.table("~/genomedk/N140_Targeting/specs/umiseq_paper/data/IMPROVEptList", header = TRUE)
prior <- readRDS("~/genomedk/N140_Targeting/specs/specs_analysis/sw_input_files/180903_prior.RDS")[, 1:4]
pon38 <- readRDS("~/genomedk/N140_Targeting/specs/umiseq_paper/reference/201020_hg38-novaseq-xgen-sporacrc-pon.RDS")
pon19 <- readRDS("~/genomedk/N140_Targeting/specs/umiseq_paper/reference/200419_novaseq-xgen-sporacrc-pon.RDS")

#Add sitemut_hg38 to plasma table (overwrite plasma)
plasma <- dplyr::left_join(plasma, wes[, c("pt_id", "sitemut_hg38")])
plasma <- plasma[!is.na(plasma$sitemut_hg38), ] #Exclude samples with no sitemut


#Add the index for the 18094*4 panel matrix
sitemut <- t(apply(pon38$coordinates, 1, function(x){
  paste( paste0(trimws(x[1]), ":", trimws(x[2]), "_"),
         paste(x[3] ,c("A", "T", "C", "G"), sep = "/"),
         sep = "")}))

plasma$i <- sapply(as.character(plasma$sitemut_hg38), function(s) which(sitemut %in% s))

#Get the pileups
pileups <- list.files("~/genomedk/IMPROVE/sporacrc/N227", recursive = T, full.names = T, pattern = "consensus.bait.pileup$")
pileups <- unique(pileups[sapply(plasma$library_id, function(l)which(grepl(l, pileups)))]) #Only get pileups needed

#Make counts assuming hg19 (otherwise, split for hg19 and hg38)
tmp <- piles_to_counts(pileups, pon19$regions)
counts <- tmp[,,1:4] + tmp[,,6:9]
mafs <- abind::abind( lapply(1:dim(counts)[[1]], function(s)(counts[s,,]/rowSums(counts[s,,]))), along = 0)
vars <- apply(counts, 2:3, var) # variance of the counts of each position based on PON
vars[ vars == 0 ] <- 1e-6


#Calculate f1 per sitemut (row in plasma)
f1 <- function(maf, pon_var, prior)prior * maf/pon_var # VAF divided by PON variance & weighted by Cosmic

plasma$score <-
  sapply(1:nrow(plasma), function(r) {
    maf0 <- mafs[ which(grepl(plasma$library_id[r], pileups)), , ][ plasma$i[r] ]  
    pon_var0 <- vars[ plasma$i[r] ]
    prior0 <- prior[ plasma$i[r] ]
    f1(maf0, pon_var0, prior0 )
  })

#Then sum score for each library_id (overwrite plasma)
plasma <-
  dplyr::group_by(plasma, library_id) %>%
  dplyr::mutate(sum = sum(score))