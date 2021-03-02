#counts array, genome version taken cared of
cruk <- readRDS("~/genomedk/matovanalysis/umiseq_analysis/R/cruk-counts.RDS") # 
dimnames(cruk)[[1]] #The pileup names are the first dim

#sample and patient info incl the pileup names
info <- read.table("~/genomedk/matovanalysis/umiseq_analysis/R/cruk-plasma-info.lst", header = T, stringsAsFactors = F)
names(info)
head(info)
table(info$sample_type)
str(info)

#Find the index in info corresponding to a slice in the cruk array.
info$i <- sapply(info$pileup, grep, dimnames(cruk)[[1]])

#For instance, sample in row 90 is slice number 27 in cruk
info[90, c("pileup", "i")] #i is 27
dimnames(cruk)[[1]][27] #The same sample

#To find cancer, we define it as sample_type with the string "CRC"
i <- unique(info[ grepl("CRC", info$sample_type), "i" ])

# To then get this (i) subset in cruk counts
cruk_cancer <- cruk[i, , ]
dim(cruk_cancer)
saveRDS(cruk_cancer, "cruk_cancer.RDS")
#Etc for other subsets

length(unique(info[ grepl("adenom", info$sample_type), "i" ]))
length(unique(info[ grepl("control", info$sample_type), "i" ]))

#How to discard INDELs from sitemuts
info$sitemut_hg38 #All the sitemuts including INDELS
p <- ".+:[[:digit:]]+_.(.)*/.(.)*" #This captures ref and alt that are longer than 1 (i.e. DEL and INS)

#So all INDELS will have 1 or more chararters in at least on of the capture groups "\\1" and "\\2"
nchar(sub(p, "\\1\\2", c("chr5:112839839_T/G", "chr5:112839838_GT/G", "chr5:112839838_G/GTA")))>0 #2 (DEL) and 3 (INS) identified

#These are SNPs (SNVs) - characters in "\\1\\2" are 0
tmp <- info[ nchar(sub(p, "\\1\\2", info$sitemut_hg38)) == 0, ]

#Remove those with no mutations (controls)
cancer_SNPs <- tmp[!is.na(tmp$sitemut_hg38), ]
#info$sitemut_hg38
#If you need the genome version find it in the pipeline parameter file param.json
#If /path/to/sample/output/pileup then /path/to/sample/param.json
#If want to access param.json locally you need change the paths - I need to change
#"/faststorage/project" to "~/genomedk/gdk-projects" like this:
tmp <- sub("/faststorage/project", "~/genomedk/gdk-projects", info$pileup)
all(file.exists(tmp)) #check that you can find the file (should be TRUE)

#The make the param paths (if /path/to/sample/output/pileup then /path/to/sample/param.json)
params <- file.path(sub("output.+", "", tmp), "param.json")
all(file.exists(params)) #Check!

#we can make a new column with genome vesion by looking for hg19 and hg38 in the param.json
#The first way beloq requires the jsonlite library and is safest.
#The second also works way also works and can be done in base R.
info$version <- sapply(params, function(p) gsub(".+((hg19)|(hg38)).+", "\\1", tolower(jsonlite::read_json(p)$reference$reference)), USE.NAMES = F)
info$version <- sapply(params, function(x)sub(".+((hg19)|(hg38)).+", "\\1", paste(readLines(x), collapse = "")), USE.NAMES = F)

cruk <- readRDS("~/projects/pileup/specs/data/cruk-counts.RDS")
info <- read.table("~/projects/pileup/specs/data/cruk-plasma-info.lst",
                   header = T, stringsAsFactors = F)

#Add the index in cruk for row in info
info$i <- sapply(info$pileup, grep, dimnames(cruk)[[1]])

#Intitate with FALSE cruk indexing array same dim as cruk except 4 bases
cruki <-
  array(rep(F, dim(cruk)[1]*dim(cruk)[2]*4),
        dim = list(dim(cruk)[1], dim(cruk)[2], 4),
        dimnames = list(dimnames(cruk)[[1]], dimnames(cruk)[[2]], dimnames(cruk)[[3]][1:4])
  )
dim(cruki)
dimnames(cruki)

#Sitemut panel in the "chr:pos_ref/alt" format
pon <- readRDS("~/genomedk/IMPROVE/call/references/201217_hg38-novaseq-xgen-sporacrc-pon.RDS")

sitemut_panel <-
  t(apply(pon$coordinates, 1, function(x){
    paste( paste0(trimws(x[1]), ":", trimws(x[2]), "_",
                  trimws(x[3]), "/", dimnames(cruki)[[3]]))}))

#The index of, say, info$sitemut_hg38[100], in every 18094*4 matrix
#(be that maf, counts, posterior or noise - if made this way) is
i0 <- grep(info$sitemut_hg38[100], sitemut_panel)

#Which is found in cruki (for patient defined by column i) and set to 1 as:
cruki[ info$i[100], , ][i0] <- 1

#Do this for all rows in info using apply or loop (as here)
p <- ".+:[[:digit:]]+_.(.)*/.(.)*" #Pattern to grap indels

for(row in 1:nrow(info)) {
  if(! is.na(info$sitemut_hg38[row]) ) {
    if(nchar(sub(p, "\\1\\2", info$sitemut_hg38[row])) == 0) {
      i0 <- grep(info$sitemut_hg38[row], sitemut_panel)
      cruki[ info$i[row], , ][i0] <- 1
    }
  }
}

sum(cruki) #114 mutations placed (INDELS removed)

#But lets manually sanity check for pt 7007 who has 2 simpel SNPs and a single DEL
info[info$pt_id %in% "7007", c("pt_id", "sitemut_hg38", "i")] #slice 56
i0 <- grep(c("chr5:112838934_C/T|chr17:7674221_G/C"), sitemut_panel)
cruki[56, , ][i0]
sum(cruki[56, , ]) == 2

#Look at the first mutation in its row/col surroundings
i0 <- which(sitemut_panel == c("chr5:112838934_C/T"), arr.ind = T)
cruki[56, , ][(i0[1]-5):(i0[1]+5), ]

saveRDS(cruki, "~/tmp/cruki.RDS")