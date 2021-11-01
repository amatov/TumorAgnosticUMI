### Tumor informed Malahnobis
# Explorative ----
library(dplyr)
#source("~/genomedk/projects/umiseq/development/tools.R") #The very latest tools
source("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/specs_analysis/sw_input_files/tools.R")

setwd("~/projects/pileup/specs/analyses/tumorinformed")

# Include cohorts and info, calculate counts, mafs
qia <- readRDS("~/projects/pileup/specs/data/qiagen-counts.RDS")
imp <- readRDS("~/projects/pileup/specs/data/improve-counts.RDS")
#pon <- readRDS("~/genomedk/IMPROVE/call/references/201217_hg38-novaseq-xgen-sporacrc-pon.RDS")
pon <- readRDS("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/umiseq_paper/reference/201217_hg38-novaseq-xgen-sporacrc-pon.RDS") # 46
cruk <- readRDS("~/projects/pileup/specs/data/cruk-counts.RDS")
info <- read.table("~/projects/pileup/specs/data/cruk-plasma-info.lst",
                   header = T, stringsAsFactors = F)
stage <- read.table("~/projects/pileup/specs/data/plasma-stage-time.lst", header = T, stringsAsFactors = F)

#Excluded
flag <- readRDS(file ="~/projects/pileup/specs/analyses/noise/2021-01-05_flagged-positions.RDS")
xx <- union_all(flag$xi, flag$xr, flag$xb) #IDSNP, REF, BALCKLIST (no need for DEL falg$xd)


#The cruk ATCG-counts and -mafs (exdlude dels)
counts <- cruk[,,1:4] + cruk[,,6:9]
mafs <- abind::abind( lapply(1:dim(counts)[[1]], function(s)(counts[s,,]/rowSums(counts[s,,]))), along = 0)
dim(mafs) #[1]   188 (192) 18094     4

qia_counts <- qia[,,1:4] + qia[,,6:9]
qia_mafs <- abind::abind( lapply(1:dim(qia_counts)[[1]], function(s)(qia_counts[s,,]/rowSums(qia_counts[s,,]))), along = 0)
dim(qia_mafs) #[1]   24 18094     4

imp_counts <- imp[,,1:4] + imp[,,6:9]
imp_mafs <- abind::abind( lapply(1:dim(imp_counts)[[1]], function(s)(imp_counts[s,,]/rowSums(imp_counts[s,,]))), along = 0)
dim(imp_mafs) #[1]   358 (378) 18094     4

#The pon ATCG variance, sd and maf
pon_counts <- pon$pon[,,1:4] + pon$pon[,,6:9] #counts
pon_mafs <- abind::abind(lapply(1:dim(pon_counts)[[1]], function(s)pon_counts[s,,]/rowSums(pon_counts[s,,], na.rm = T)), along = 0)

var <- apply(pon_mafs, 2:3, var, na.rm = T) #variance of mafs per substitution across the pons (46! for now)
sum( var == 0) #13459 including xx's
hist(log10(var))
var[var == 0] <- 1E-11 #Seem ok
h <- hist(log10(var[-xx]), plot = F)
col <- c("red", rep("blue", length(h$breaks)-1))
plot(h, col = col, main = "PON core substitution mafs", xlab = "Log10(var)")

stdev <- apply(pon_mafs, 2:3, sd, na.rm = T) #sd per substitution across the pons (46! for now)
hist(log10(stdev))
stdev[ stdev == 0] <- 1E-6 #10x lower than min
h <- hist(log10(stdev[-xx]), plot = F)
col <- c("red", rep("blue", length(h$breaks)-1))
plot(h, col = col, main = "PON substitution mafs (core)", xlab = "Log10(stdev)")


#Calculate mahal array (as below) for all counts
#mahal_old <- function(maf, weight, prior = 1, e = 0)prior * maf/(weight + e) # VAF divided by PON 1/weight and prior
#Generalized score per substition i maf_i*w1_i*w2_i, ... where, say the weights 1/std and prior is supplied in ...:
#mahal(maf, 1/std, prior)
mahal <- function(maf, ...)Reduce(`*`, list(maf, ...))

#score_old <- abind::abind( lapply(1:dim(mafs)[[1]], function(s)mahal_old(mafs[s,,], stdev, e = 0)), along = 0) #W stdev
#range(score_old, na.rm = T)
score <- abind::abind( lapply(1:dim(mafs)[[1]], function(s)mahal(mafs[s,,], 1/stdev)), along = 0) #W stdev
range(score, na.rm = T)
#saveRDS(score, "~/tmp/mahal-score-cruk.RDS")
#The other scores
pon_score <- abind::abind( lapply(1:dim(pon_mafs)[[1]], function(s)mahal(pon_mafs[s,,], 1/stdev)), along = 0)
qia_score <- abind::abind( lapply(1:dim(qia_mafs)[[1]], function(s)mahal(qia_mafs[s,,], 1/stdev)), along = 0)
imp_score <- abind::abind( lapply(1:dim(imp_mafs)[[1]], function(s)mahal(imp_mafs[s,,], 1/stdev)), along = 0)

hist(log10(pon_score[45,,][-xx] + 1E-4 ), xlim = c(-4, 4)) #Add 1-E4 to display the 0s
#Notice a few with signal especially in the range 1 - 20-ish
hist(log10(qia_score[1,,][-xx] + 1E-4 ), xlim = c(-4, 4), main = "Single scores example (control QIA1, core)", xlab = "Log10(score + 1E-4)")
hist(log10(qia_score[1,,][-xx]), breaks = 100, xlim = c(-0.1, 1.5), ylim = c(0, 250), main = "Single scores example (control QIA1, core)", xlab = "Log10(score)")
sum(qia_score[1,,][-xx] > 0.1)
hist(log10(imp_score[10,,][-xx] + 1E-4 ), xlim = c(-4, 4))

#Make image hitmaps of multiple samples
source("~/projects/utility_functions/image_plot.R")
#Cut scores into the intevals
cuts <- c(-Inf, 1E-4, 1E-2, 1, 5, 11, 16, 50, +Inf) #Good for pon mafs based dispersion
# cuts <- c(-Inf, 1E-5, 1E-4, 1E-3, 0.01, 0.1, 1, +Inf) #Good for pon counts based dispersion
#For each sample count observations in cuts
pon_t <- apply(pon_score, 1, function(x) table(cut(x[-xx], breaks = cuts)))
cruk_t <- apply(score, 1, function(x) table(cut(x[-xx], breaks = cuts)))
qia_t <- apply(qia_score, 1, function(x) table(cut(x[-xx], breaks = cuts)))
imp_t <- apply(imp_score, 1, function(x) table(cut(x[-xx], breaks = cuts)))


bp <- c(0, 5, 20, 50, 200, 1000, 5000, 20000, 30000, 40000, 100000)
dev.off()
pon_t
image_plot(t(as.matrix(pon_t)), legend = T, xlab = rownames(pon_t), col = terrain.colors(11), breaks = bp, legend_title = "Freq", plot_title = "PON single scores")
image_plot(t(as.matrix(cruk_t)), legend = T, xlab = rownames(cruk_t), col = terrain.colors(10), breaks = bp, , legend_title = "Freq", plot_title = "CRUK single scores")
png("figures/single-score-w1-qia.png", res = 450, width = 3000, height = 2200)
image_plot(t(as.matrix(qia_t)), legend = T, xlab = rownames(qia_t), col = terrain.colors(10), breaks = bp, legend_title = "Freq", plot_title = "QIA controls single scores")
dev.off()
image_plot(t(as.matrix(imp_t)), legend = T, xlab = rownames(imp_t), col = terrain.colors(10), breaks = bp, legend_title = "Freq", plot_title = "IMP single scores")

#Separate imp_t into pre and post (most should be negative)
imp_prop <- stage$library_id[stage$op_time_cat == "-1"]
imp_poop <- stage$library_id[stage$op_time_cat %in% c("2", "30")]
imp_prop_t <- imp_t[, unlist(sapply(imp_prop, grep, x = dimnames(imp)[[1]]))]
imp_poop_t <- imp_t[, unlist(sapply(imp_poop, grep, x = dimnames(imp)[[1]]))]
image_plot(t(as.matrix(imp_prop_t)), legend = T, xlab = rownames(imp_t), col = terrain.colors(10), breaks = bp, legend_title = "Freq", plot_title = "IMP preOP single scores")
image_plot(t(as.matrix(imp_poop_t)), legend = T, xlab = rownames(imp_t), col = terrain.colors(10), breaks = bp, legend_title = "Freq", plot_title = "IMP postOP single scores")

#Also Cruk
cruk_prop <- unique(info$library_id[info$sample_type %in% "CRC pre-OP"])
cruk_ctrl <- unique(info$library_id[grepl("control", info$sample_type)])
cruk_prop_t <- cruk_t[, unlist(sapply(cruk_prop, grep, x = dimnames(counts)[[1]])) ] 
cruk_ctrl_t <- cruk_t[, unlist(sapply(cruk_ctrl, grep, x = dimnames(counts)[[1]])) ]

png("figures/single-score-w1-cruk-prop.png", res = 450, width = 3000, height = 2100)
image_plot(t(as.matrix(cruk_prop_t)), legend = T, xlab = rownames(imp_t), col = terrain.colors(10), breaks = bp, legend_title = "Freq", plot_title = "CRUK preOP single scores")
dev.off()
png("figures/single-score-w1-cruk-ctrl.png", res = 450, width = 3000, height = 2100)
image_plot(t(as.matrix(cruk_ctrl_t)), legend = T, xlab = rownames(imp_t), col = terrain.colors(10), breaks = bp, legend_title = "Freq", plot_title = "CRUK ctrl single scores")
dev.off()




#Add the index in cruk for row in info for easy lookup
info$i <- sapply(info$pileup, grep, dimnames(cruk)[[1]])

#Sitemut panel in the "chr:pos_ref/alt" format
sitemut_panel <-
  t(apply(pon$coordinates, 1, function(x){
    paste( paste0(trimws(x[1]), ":", trimws(x[2]), "_",
                  trimws(x[3]), "/", dimnames(counts)[[3]]))}))
dim(sitemut_panel) #[1] 18094     4

#Add panel matrix index of mutations (INDELs and NAs get NA)
info$j <- sapply(info$sitemut_hg38, function(x) ifelse(!is.na(x), grep(pattern = x, x = sitemut_panel), NA))
# info[ is.na(info$j), "sitemut_hg38" ]
# #Intitate with 0 cruk indexing array same dim as cruk counts bases A T C G
# cruki <-
#   array(rep(0, dim(counts)[1]*dim(counts)[2]*dim(counts)[3]),
#         dim = dim(counts),
#         dimnames = dimnames(counts)
#   )
# 
# #Do this for all rows in info using apply or loop (as here)
# # p <- ".+:[[:digit:]]+_.(.)*/.(.)*" #Pattern to grap indels (removed below)
# p <- ".+_./.$" #Pattern to grap SNVs
# 
# for(row in 1:nrow(info)) {
#   if(grepl(p, info$sitemut_hg38[row])) {
#     i0 <- grep(info$sitemut_hg38[row], sitemut_panel)
#     cruki[ info$i[row], , ][i0] <- 1
#   }
# }
# sum(cruki) #168 mutations placed (INDELS removed)
#
# #But lets manually sanity check for pt 7007 who has 2 simpel SNPs and a single DEL
# info[info$pt_id %in% "7007", c("pt_id", "sitemut_hg38", "i")] #slice 56
# i0 <- grep(c("chr5:112838934_C/T|chr17:7674221_G/C"), sitemut_panel)
# cruki[56, , ][i0]
# sum(cruki[56, , ]) == 2
# 
# #Look at the first mutation in its row/col surroundings
# i0 <- which(sitemut_panel == c("chr5:112838934_C/T"), arr.ind = T)
# cruki[56, , ][(i0[1]-5):(i0[1]+5), ]
# saveRDS(cruki, "~/tmp/cruki.RDS")

# #What are the MAFS and scores for the preOP with SNV in cruk
# #These CRUK indeces are preOP cancers
# i <- 
#   info %>%
#   filter(sample_type %in% c("CRC high ctDNA", "CRC pre-OP"), 
#          grepl(".+_./.$", sitemut_hg38)) %>% 
#   .$i %>% 
#   unique()
# 
# res <- c()
# 
# for(x in 1:length(i)) {
#   sitemut0 <- grep(".+_./.$", info$sitemut_hg38[info$i == i[x]], value = T)
#   library_id0 <- rep(info[info$i == i[x], "library_id"][1], length(sitemut0))
#   maf0 <- mafs[i[x],,] [ cruki[i[x],,] == 1 ]
#   score0 <- score[j[x],,] [ cruki[i[x],,] == 1 ]
#   res <- rbind(res, cbind(library_id0, sitemut0, maf0, score0 ) )
# }
# score_cruk <- as.data.frame(res, stringsAsFactors = F)
# plot(as.numeric(score_cruk[,"score0"])~as.numeric(score_cruk[,"maf0"]), log = "y")
# 
# #Summerize per sample
# pdata <- 
#   score_curk %>% 
#   mutate(maf = as.numeric(maf0),
#          score = as.numeric(score0)) %>% 
#   select(-score0, -maf0) %>% 
#   group_by(library_id0) %>% 
#   summarise(maf = mean(maf), 
#             score = sum(score)) %>% 
#   as.data.frame()
# filter(pdata, score==0) %>% .$library_id0 %>%  paste(collapse = ",")
# par(mfrow = c(1,2))
# with(pdata, plot(score ~ maf, log = "y", main = "", xlab = "mean_maf"))
# with(pdata, hist(log(score), main = ""))
# mtext(side = 3, line = -2, "Plasma sum-of-maf/var scores\n(68 preOP samples - of which 11 are null)", outer = T)
# dev.off()



#Lets go through preOP samples and get the score for sample i vs other cancer and controls
# The cancer patients (take those with )
ptids <- unique(info[info$sample_type %in% "CRC pre-OP" & grepl("chr", info$sitemut_hg38), "pt_id"])
ica <- unique(info$i[info$pt_id %in% ptids]) #All preOP ids
ino <- unique(info$i[grepl("(adenoma)|(control)", info$sample_type)]) #All control ids

score_res <- list()
for(pt_id in ptids) {
  
  sitemut0 <- info$sitemut_hg38[info$pt_id %in% pt_id] #Sitemuts of pt_id
  imut0 <- which(sitemut_panel %in% sitemut0)#Index of mutations for sample i0 in the 18094*4 matrix ignoring INDELS 
  i0 <- info$i[info$pt_id %in% pt_id][1] #The panel index  - take the first as they are identical
 
  
  score_sample0 <- sum(mahal(mafs[i0,,][imut0], vars[imut0]), na.rm = T)
  score_cancer0 <- sapply(ica, function(x) sum(mahal(mafs[x,,][imut0], vars[imut0]), na.rm = T))
  score_control0 <- sapply(ino, function(x) sum(mahal(mafs[x,,][imut0], vars[imut0]), na.rm = T))
  score_pon0 <- sapply(1:dim(pon_mafs)[1], function(x) sum(mahal(pon_mafs[x,,][imut0], vars[imut0]), na.rm = T))
  

  boots <- replicate(1000, sample(1:length(vars), length(imut0)), simplify = F)
  score_boots0 <- sapply(boots, function(x) sum(mahal(mafs[i0,,][x], vars[x]), na.rm = T))
    
  score_res[[as.character(pt_id)]] <-
    list(sample = score_sample0,
         cancer = score_cancer0,
         control = score_control0,
         pon = score_pon0,
         boot = score_boots0,
         sitemuts = imut0)
  
}
par(mar=c(6,4,2,2))
plot.new()
plot.window(ylim = c(1E-4, max(unlist(score_res))), xlim = c(0, length(ptids)), log = "y")
axis(2)
x=1
for(s in names(score_res)) {
  
  points(x = rep(x-0.25, length(score_res[[s]]$control)), y = score_res[[s]]$control + 1E-4)
  points(x = rep(x-0.25, length(score_res[[s]]$pon)), y = score_res[[s]]$pon + 1E-4, pch = 0)
  points(x = rep(x+0.25, length(score_res[[s]]$cancer)), y = score_res[[s]]$cancer + 1E-4, col = "red") 
  points(x = x, y = score_res[[s]]$sample + 1E-4 , col = "blue", cex = 2)
  p <- sprintf("%.3f",(sum(score_res[[s]]$boot > score_res[[s]]$sample)/1000))
  text(p, x = x, y = 1E-5, srt = 90, cex = 0.75, xpd = T)
  x <- x + 1
}

#TODO substract mean score in normals before p value estimation



#Plot differences
plot.new()

plot.window(ylim = c(1E-4, max(unlist(score_res))), xlim = c(0, length(ptids)), log = "y")
axis(2)
x=1
for(s in names(score_res)) {
  
  points(x = rep(x-0.25, length(score_res[[s]]$control)), y = score_res[[s]]$control + 1E-4)
  points(x = rep(x+0.25, length(score_res[[s]]$cancer)), y = score_res[[s]]$cancer + 1E-4, col = "red") 
  points(x = x, y = score_res[[s]]$sample + 1E-4 , col = "blue", cex = 2)
  
  x <- x + 1
}


# Make sw-AND posteriors and sw-AND baysian factors for cruk ----
# TODO: Redo with bamclipping - see .../project/test/bamclip script
# source("~/genomedk/IMPROVE/sporacrc/current_version/tools.R")
# 
# pileups <- sub("/faststorage/project/", "/home/mhra/genomedk/gdk-projects/", dimnames(cruk)[[1]])
# 
# #Get the param.json files using pileups as scaffold
# params <- file.path(sub("/output(?!(.*output.*))", "", dirname(pileups), perl = T), "param.json") #Assume param is in "the last" output
# all(file.exists(params))
# 
# #Extract reference version from param file (field reference$reference)
# version <- sapply(params, function(p) gsub(".*(hg19|hg38).*", "\\1", tolower(jsonlite::read_json(p)$reference$reference)), USE.NAMES = F)
# 
# 
# #Use pileups_to_counts to get counts (sw_piles to get both posterior and counts) - takes a while with 100s samples
# pon19 <- readRDS("~/genomedk/N140_Targeting/specs/umiseq_paper/reference/200419_novaseq-xgen-sporacrc-pon.RDS")
# pon38 <- readRDS("~/genomedk/N140_Targeting/specs/umiseq_paper/reference/201217_hg38-novaseq-xgen-sporacrc-pon.RDS")
# 
# #TODO: get both posterior and bf out
# 
# for(v in c("hg19", "hg38")) {
#   if(v == "hg19") {
#     post19 <-
#       lapply(pileups[version %in% v], function(p) {
#         sw_piles(pileup = p,
#                  pon = pon19$pon,
#                  regions = pon19$regions,
#                  prior = 0.5, model = "AND", rho = NULL)$posterior
#       })
# 
#     names(post19) <- gsub(".*(N[[:digit:]]{3}-[[:digit:]]+)_.*", "\\1", pileups[version %in% v])
#   }
#   if(v == "hg38") {
#     post38 <-
#       lapply(pileups[version %in% v], function(p) {
#         sw_piles(pileup = p,
#                  pon = pon38$pon,
#                  regions = pon38$regions,
#                  prior = 0.5, model = "AND", rho = NULL)$posterior
#       })
# 
#     names(post38) <- gsub(".*(N[[:digit:]]{3}-[[:digit:]]+)_.*", "\\1", pileups[version %in% v])
#   }
# }
# 
# str(post)
# names(post)
# 
# post <- c(post19, post38)
# post <- abind::abind(post, along = 0)
# saveRDS(post, "cruk-AND.RDS")

post <- readRDS("~/projects/pileup/specs/data/cruk-AND.RDS")


#Use cruki above to extract sitemut posteriors similar as we extract mahal sum
# What are the posteriors preOP of SNV in cruk - Note ordering in post differ from score
# Slice index of WES mutations in preOP

#Add the index in post for row in info for easy lookup using library_id
info$posti <- sapply(info$library_id, function(x)grep(paste0("^", x, "$"), dimnames(post)[[1]]))
info[,c("i", "j", "posti")]
j <- 
  info %>%
  filter(sample_type == "CRC pre-OP", 
         grepl(".+_./.$", sitemut_hg38)) %>% 
  .$posti %>% 
  unique()
info[ c("i", "posti")]
res <- c()
for(x in 1:length(j)) {
  sitemut0 <- grep(".+_./.$", info$sitemut_hg38[info$i == j[x]], value = T)
  library_id0 <- rep(info[info$i == j[x], "library_id"][1], length(sitemut0))
  
  maf0 <- mafs[ j[x],,] [ cruki[j[x],,] == 1 ]
  count0 <- counts[ j[x],,] [ cruki[j[x],,] == 1 ]
  post0 <- post[[ j[x] ]] [ cruki[j[x],,] == 1 ]
  
  res <- rbind(res, cbind(library_id0, sitemut0, maf0, count0, post0 ) )
}
#TODO for score and post: go through info sample instead and parse per lib and sitemut. 
#Now there are issue with recurrent sitemuts. 
# 
filter(info, sitemut_hg38 == "chr5:112839726_C/T")

post_cruk <- as.data.frame(res, stringsAsFactors = F)

full_join(score_cruk, post_cruk, by = c("library_id0", "sitemut0"))

plot(as.numeric(score_cruk[,"score0"])~as.numeric(score_cruk[,"maf0"]), log = "y")


# Make sw and scores ----
# Scores an posteriors from above

#Initiate result df
res <- info[!is.na(info$j), c("library_id", "sitemut_hg38", "i", "j")]


## The cruk index slice i can be used as posterior index also - they are identical
#all(res$i == sapply(paste0("^", res$library_id, "$"), grep, dimnames(post)[[1]]))

#For each library_id lookup mahal and post and place into res0
res0 <- c()
for(row in 1:nrow(res)) {
  i0 <- res$i[row] #The data array slice index of sample
  j0 <- res$j[row] #The panel index of sitemut or NA if not SNV
  res0 <- rbind(res0, c(score[i0,,][j0], post[i0,,][j0], mafs[i0,,][j0], counts[i0,,][j0]))
}
res <- cbind(res, res0)  
names(res) <- c("library_id", "sitemut_hg38", "i", "j", "mahal", "sw", "maf", "count")  
res <- 
  res %>%
  dplyr::group_by(library_id) %>% 
  dplyr::mutate(score_mahal = sum(mahal),
                score_minpost = min(sw),
                score_fishpost = -2*sum(log(sw)))
dev.off()

source("~/projects/utility_functions/recoder.R")
pdata <- 
  dplyr::select(res, library_id, starts_with("score") ) %>%
  dplyr::mutate(n = n()) %>% 
  dplyr::distinct()

pchs <- as.character(pdata$n) 
cols <- recoder(pdata$score_minpost < 0.05, match_list = list(0:1, c("black", "red")))
plot(log10(pdata$score_mahal + 1E-5), log10(pdata$score_fishpost), col = cols, pch = pchs)
filter(pdata, score_mahal == 0)
arrange(res, maf) %>% as.data.frame()

#Make som empirical P values for samples and control
f1 <- function(v, data, n, N, fun) {
  #Estimate empirical p-value of number v by applying function fun on n draws in data
  #using N boots. When v = NULL, return distribution, else (r+1)/(N+1) - r fraction of better boots than v
  boots <- replicate(N, fun(sample(data, n)))
  if(is.null(v))return(boots)
  (sum(boots > v)+1) / (N + 1) #Fraction of boots with greater teststat/score
}
#For bootstrapping we need to exclude all flagged (not dels not considered)
flag <- readRDS(file ="~/projects/pileup/specs/analyses/noise/2021-01-05_flagged-positions.RDS")
xx <- dplyr::union_all(flag$xi, flag$xr, flag$xb)

#Takes some minutes for 100+ samples and N=100000 boots
#
# TODO - CHECK Seems like we include the excluded, not include (xx should be -xx ??)
pboots <- 
  dplyr::select(res, library_id, i, starts_with("score") ) %>%
  dplyr::mutate(n = n()) %>% 
  dplyr::distinct() %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(p_mahal = f1(score_mahal, score[i,,][xx], n, 100000, fun = function(x)sum(x, na.rm = T) ), 
                p_fishpost = f1(score_fishpost, post[i,,][xx], n, 100000, fun = function(x)-2*sum(log(x), na.rm = T)))


cols <- recoder(pboots$score_minpost < 0.05, match_list = list(0:1, c("black", "red")))
plot(log10(pboots$p_mahal) ~ log10(pboots$p_fishpost), pch = as.character(pboots$n), col = cols)
abline(0,1, col = "red")
abline(v = log10(0.05), col = "blue")
abline(h = log10(0.05), col = "blue")

set.seed(8)
f1(0.0108, score[125,,][xx], 3, 100000, function(x)sum(x, na.rm = T) )
f1(136, post[125,,][xx], 3, 100000, function(x)-2*sum(log(x), na.rm = T) )

p <- sprintf("%.3f",(sum(score_res[[s]]$boot > score_res[[s]]$sample)/1000))



# Testing  ----------------------------------------------------------------
# Take the best from above to generate test stats (e.g. ROC AUC) for:
# 1) SW_best 2) SW_comb 3) mahal based on cruk

library(dplyr)
source("~/genomedk/projects/umiseq/development/tools.R") #The very latest tools
source("~/projects/utility_functions/auc.R") #use auc_mw()
source("~/projects/utility_functions/recoder.R")
setwd("~/projects/pileup/specs/analyses/tumorinformed")

# Include cohorts and info, calculate counts, mafs
cruk <- readRDS("~/projects/pileup/specs/data/cruk-counts.RDS")
pon <- readRDS("~/genomedk/IMPROVE/call/references/201217_hg38-novaseq-xgen-sporacrc-pon.RDS")
info <- read.table("~/projects/pileup/specs/data/cruk-plasma-info.lst", header = T, stringsAsFactors = F)

table(info$sample_type)
info %>% 
  filter(is.na(sitemut_hg38), !grepl("control", sample_type))
         


#Excluded
flag <- readRDS(file ="~/projects/pileup/specs/analyses/noise/2021-01-05_flagged-positions.RDS")
xx <- union(flag$xi, union(flag$xr, flag$xb)) #IDSNP, REF, BALCKLIST (no need for DEL flag$xd)


#The cruk ATCG-counts and -mafs (exclude dels in cols 5 and 10)
counts <- cruk[,,1:4] + cruk[,,6:9]
mafs <- abind::abind( lapply(1:dim(counts)[[1]], function(s)(counts[s,,]/rowSums(counts[s,,]))), along = 0)
dim(mafs) #[1]   201 18094     4

#The pon ATCG variance, sd and maf
pon_counts <- pon$pon[,,1:4] + pon$pon[,,6:9] #counts
pon_mafs <- abind::abind(lapply(1:dim(pon_counts)[[1]], function(s)pon_counts[s,,]/rowSums(pon_counts[s,,], na.rm = T)), along = 0)
stdev <- apply(pon_mafs, 2:3, sd, na.rm = T) #sd per substitution across the pons (46! for now)
stdev[ stdev == 0] <- 1E-6 #10x lower than min (see previously)

#Calculate mahal array (as below) for all counts
mahal <- function(maf, ...)Reduce(`*`, list(maf, ...))

score <- abind::abind( lapply(1:dim(mafs)[[1]], function(s)mahal(mafs[s,,], 1/stdev)), along = 0) #W stdev

#AND posteriors 
post <- readRDS("~/projects/pileup/specs/data/cruk-AND.RDS")

#Sitemut panel in the "chr:pos_ref/alt" format
sitemut_panel <-
  t(apply(pon$coordinates, 1, function(x){
    paste( paste0(trimws(x[1]), ":", trimws(x[2]), "_",
                  trimws(x[3]), "/", dimnames(counts)[[3]]))}))
dim(sitemut_panel) #[1] 18094     4

#Add panel matrix index of mutations (INDELs and NAs get NA and is remove below)
info$sitemuti <- sapply(info$sitemut_hg38, function(x) ifelse(!is.na(x), grep(pattern = x, x = sitemut_panel), NA))

#Add the cruk array index: cruki
info$cruki <- sapply(info$pileup, grep, dimnames(cruk)[[1]])

#Add the array index for the SW posteriors (labelled by library_id e.g. N289-124)
posti <- sapply(info$library_id, function(x)grep(paste0("^", x, "$"), dimnames(post)[[1]]))
posti[sapply(posti, length) == 0] <- NA #Some info samples were not in post so replace with NA
info$posti <- unlist(posti)

#Finally filter the minimal info for complete control and case data - other info 
#e.g. comorbidities can be include if needed - here just distinguish cancer/control.
cases <-
  info %>% 
  dplyr::filter(!is.na(sitemuti), #panel index for mutation (i.e. is SNV)
                !is.na(cruki), #array index in cruk data
                !is.na(posti), #array index in sw posterior data
                !sitemuti %in% xi) %>%   #remove a few flagged (idsnp regions)
  dplyr::mutate(cancer = 2) %>% #recode - set adenomas to 1 and controls to 0
  dplyr::select(library_id, cancer,  sitemuti, cruki, posti)

controls <- 
  info %>% 
  dplyr::filter(grepl("adenom|control", sample_type), #Defines the controls (put adenomas in for now)
                !is.na(cruki), #array index in cruk data
                !is.na(posti)) %>%  #array index in sw posterior data
  dplyr::mutate(cancer = ifelse(grepl("adenom", sample_type), 1, 0)) %>% 
  dplyr::select(library_id, cancer, sitemuti, cruki, posti)

unique(cases$library_id)
table(controls$cancer)

#Find ROC AUC by two measures: 1) looking up the sitemut sets known, 2) sampling random,
#(3) sampling by COSMIC prior ?)
# The test statistic can be the score (mahal, sw_best, sw_comb) or an empirical p value

#For the cases add maf, counts and score and post per mutation
res0 <- c()
for(row in 1:nrow(cases)) {
  i0 <- cases$sitemuti[row] #The panel index of sitemut or NA if not SNV
  j0 <- cases$cruki[row] #The data array slice index of sample
  k0 <- cases$posti[row] #The posterior array index of sample
  res0 <- rbind(res0, c(score[j0,,][i0], post[k0,,][i0], mafs[j0,,][i0], counts[j0,,][i0]))
}
res0 <- data.frame(res0)
names(res0) <- c("mahal", "post", "maf", "count")
cases <- cbind(cases, res0)  

cases <- 
  cases %>%
  dplyr::group_by(library_id) %>% 
  dplyr::mutate(score_mahal = sum(mahal),
                score_minpost = min(post + 1E-30),
                score_fishpost = -2*sum(log(post + 1E-30)))


#This concludes the combined scores of the positive cases. How extreme are they ?
#We can do an empirical p value (since we dont know an exact distribtion).
#On the whole cohort case/ctrl we can calculate AUC. 
#By just looking at the same sets in controls we 

dev.off()

pdata <- 
  dplyr::select(cases, library_id, count, starts_with("score") ) %>%
  dplyr::mutate(n = n(),
                total_counts = sum(count)) %>% 
  dplyr::select(-count) %>% 
  dplyr::distinct()

pchs <- as.character(pdata$n) 
cols <- recoder(pdata$score_minpost < 0.05, match_list = list(0:1, c("black", "blue")))
cols[ pdata$total_counts < 1 ] <- "red"
cuts <- cut(pdata$total_counts, breaks = c(-Inf, 0, 2, 5, 10, 50, Inf))
cex <- recoder(cuts, list(sort(unique(cuts)), c(1, 1.2, 1.4, 1.6, 1.8, 2)))
plot(log10(pdata$score_mahal + 1E-5), log10(pdata$score_fishpost), col = cols, pch = pchs, cex = cex, 
     main = "Properties of integrated callers", ylab = "Log10(SW_comb)", xlab = "Log10(Mahal + 1E-5)")
filter(pdata, score_mahal == 0)
arrange(cases, maf) %>% as.data.frame()

#Make som empirical P values for samples and control
f1 <- function(v, data, n, N, fun) {
  #Estimate empirical p-value of number v by applying function fun on n draws in data
  #using N boots. When v = NULL, return distribution, else (r+1)/(N+1) - r fraction of better boots than v
  boots <- replicate(N, fun(sample(data, n)))
  if(is.null(v))return(boots)
  (sum(boots > v)+1) / (N + 1) #Fraction of boots with greater teststat/score
}

#Takes 5 minutes for 100+ samples and N=100000 boots
t0 <- Sys.time()
pboots <- 
  dplyr::select(cases, library_id, cruki, posti, starts_with("score") ) %>%
  dplyr::mutate(n = n()) %>% 
  dplyr::distinct() %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(p_mahal = f1(score_mahal, score[cruki,,][-xx], n, 100000, fun = function(x)sum(x, na.rm = T) ),
                p_fishpost = f1(score_fishpost, post[posti,,][-xx], n, 100000, fun = function(x)-2*sum(log(x), na.rm = T)))
t1 <- Sys.time()-t0

#Distributions of mahal and fishpost examples. Repeat to see the reproducibility.
#For a high maf sample N289-223, intermediate N289-46 and a low N289-272
sample_n(ungroup(cases), 20) # Find some representative ones
libraryids <- c("N289-223", "N289-46", "N289-144", "N289-272")
opar <- par(no.readonly = T)
par(opar)
dev.off()

for(lib in libraryids) {
  data0 <- 
    ungroup(cases) %>% 
    dplyr::filter(library_id == lib) %>% 
    dplyr::mutate(n = n(), 
                  load = mean(maf, na.rm = T)) %>% 
    dplyr::select(library_id, cruki, posti, n, load, score_mahal, score_fishpost) %>% 
    dplyr::distinct()
  
  #For this we need to remove personal SNP to not sample these (rarely) - add to xi as xx
  xx0 <- union(which( mafs[data0$cruki,,] > 0.35 ), xx)
  
  mdistr0 <- f1(v = NULL, data = score[data0$cruki,,][-xx0], n = data0$n, N = 100000, fun = function(x)sum(x, na.rm = T))
  pdistr0 <- f1(v = NULL, data = (post[data0$posti,,][-xx0])+1E-30, n = data0$n, N = 100000, fun = function(x)-2*sum(log(x), na.rm = T))
  #The original individual scores
  # plot(log10(sort(post[data0$posti,,][-xx])[1:500]), pch = 20, col = "blue")
  # plot(log10(sort(score[data0$cruki,,][-xx], decreasing = T)[1:500]), pch = 20, col = "blue")
  
  ymax <- ceiling(max(log10(sort(mdistr0, decreasing = T))[1], log10(sort(pdistr0, decreasing = T))[1] ))
  sub <- paste0("load:", round(data0$load, 4), "   muts:", data0$n, "   mahal:", round(data0$score_mahal, 4),
                "   sw_comb:", round(data0$score_fishpost, 4))
  plot(log10(sort(mdistr0, decreasing = T))[1:1000], pch = 20, col = "blue", ylim = c(0, ymax), 
       cex = 0.5, ylab = "Log10(score)", xlab = "Top 1E3 (N = 1E5)", sub = sub)
  title(main = lib, line = 0)
  points(log10(sort(pdistr0, decreasing = T))[1:1000], pch = 20, col = "red", cex = 0.5)
  legend("topright", c("mahal", "sw_comb"), col = c("blue", "red"), pch = 20)
  abline(h = c(log10(data0$score_mahal), log10(data0$score_fishpost)), col = c("blue", "red") )
  # plot(log10(sort(pdistr0, decreasing = T))[1:500], pch = 20, col = "red")
}

#Or show stability by repeating sampling for a single score (a low and a high)
set.seed(8)
sim <- 
  rbind(replicate(5, expr = {f1(3, score[125,,][-xx], 3, 100000, function(x)sum(x, na.rm = T) )}),
        replicate(5, expr = {f1(150, score[125,,][-xx], 3, 100000, function(x)sum(x, na.rm = T) )}),
        replicate(5, expr = {f1(3, post[125,,][xx], 3, 100000, function(x)-2*sum(log(x), na.rm = T))}),
        replicate(5, expr = {f1(150, post[125,,][xx], 3, 100000, function(x)-2*sum(log(x), na.rm = T))})
  )
dev.off()
plot.new()
plot.window(ylim = log10(range(sim)), xlim = c(0,5))
col = c("blue", "blue", "red", "red")
for(r in 1:nrow(sim)){

    points(y = log10(sim[r,]), jitter(rep(r, 5), 10/r), pch = 20, cex = 2, col = rep(col[r], 5))
}
axis(2)
mtext(c("mh3", "mh150", "swc3", "swc150"), side = 1, at = 1:4)
title(ylab = "Log10(p*)", main = "p* stability", sub = "(N289-245; N=1E5; n=3; core)")
#Plot empirical p-value
pchs <- as.character(pboots$n) 
cols <- recoder(pboots$score_minpost < 0.05, match_list = list(0:1, c("black", "blue"))) #p > 0.05 sw_best black
#Recode mean ctDNA load into cex size - first get mean_maf
load <-
  left_join(ungroup(pboots), cases, by = "library_id") %>% 
  dplyr::group_by(library_id) %>% 
  dplyr::summarise(load = mean(maf))
#Recode by quantiles, and add back to pboots data in correct order. 
load$cex <- recoder(load$load, list(NA, c(1, 1.5, 2)), q = c(0, 0.3, 0.66, 1))
cex <- dplyr::left_join(pboots, load) %>% .$cex
plot(y = log10(pboots$p_mahal), x = log10(pboots$p_fishpost), pch = as.character(pboots$n), 
     col = cols, cex = cex, ylim = c(-5, 0), xlim = c(-5, 0),
     main = "P* (1E5), mutations and ctDNA load", ylab = "Log10(P*_mahal)", xlab = "Log10(P*_sw_comb)")
abline(0,1, col = "red")
abline(v = log10(0.05), col = "blue")
abline(h = log10(0.05), col = "blue")


#ROC AUC
pboots
cases
controls_tmp <- controls
controls <- controls[controls$cancer == 0, ] #Remove adenomas

simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
    #Plot as points(TPR ~ 1 - FPR)
  }

#1. 
res <- 
  sapply( unique(cases$library_id), function(case) {
  i0 <- cases$sitemuti[cases$library_id == case] #The sitemut set for case
  sapply(controls$cruki, function(control) {
    sum(score[control,,][i0], na.rm = T)
  }
  )
}
)
#Get a control*cases-mutation-set matrix
dimnames(res)[[1]] <- controls$library_id
length(res) #3666

# #Let sanitty check red[2,2] "N289-198*N289-224" (res[2,2])
# i0 <- cases$sitemuti[cases$library_id=="N289-224"]
# sum(score[controls$cruki[controls$library_id == "N289-198"],,][i0]) == res[2,2] #T

#2 Same for sw_comb

res2 <- 
  sapply( unique(cases$library_id), function(case) {
    i0 <- cases$sitemuti[cases$library_id == case] #The sitemut set for case
    sapply(controls$cruki, function(control) {
      -2*sum(log(post[control,,][i0] + 1E-30), na.rm = T)
    }
    )
  }
  )
#Get a control*cases-mutation-set matrix
dimnames(res2)[[1]] <- controls$library_id
length(res) #3666



pred <- as.vector(res) #Here we still have adenomas as controls!
labels <- rep(0, length(pred))
pred <- c(pred, pboots$score_mahal) #Here we have those without signal!
labels <- c(labels, rep(1, nrow(pboots)))
auc1 <- auc_mw(labels, pred) #[1] 0.8455663
tmp <- simple_roc(labels, pred)
plot(tmp$TPR ~ (1-tmp$FPR) , type = "l", col = "blue")


pred2 <- as.vector(res2) #Here we still have adenomas as controls!
labels2 <- rep(0, length(pred2))
pred2 <- c(pred2, pboots$score_fishpost) #Here we have those without signal!
labels2 <- c(labels2, rep(1, nrow(pboots)))
auc2 <- auc_mw(labels2, pred2) #[1] 0.7977097
tmp2 <- simple_roc(labels2, pred2)
points(tmp2$TPR ~ (1-tmp2$FPR) , type = "l", col = "red" )






