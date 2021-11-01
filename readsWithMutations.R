

# raw <- read.table("C16A07111D_cfdna_N289_70_consensus.txt", sep = "\t", header = F) #Read data
# raw <- read.table("pt4mamma_cfdna_N289_295_consensus.txt", sep = "\t", header = F) #Read data
# raw <- read.table("C16A07111D_cfdna_N289_70_consensus.txt", sep = "\t", header = F) #Read data


cruk_files <- readLines("~/genomedk/matovanalysis/umiseq_analysis/R/cruk-pileups.lst") 

#Problems with some cruk bam files that are still only mapped top HG19 - we will only use the hg38s
# Example of naming that identifies hg38 and hg19
# good <- "/faststorage/project/PolyA/BACKUP/CRUK/plasma/N289/451145/S07A09782D_cfdna_N289-181___210422140851sep/output/S07A09782D_cfdna_N289-181_consensus.bait.pileup"
# bad <- "/faststorage/project/PolyA/BACKUP/CRUK/plasma/N289/6597/C30A06597D_cfdna_N289-12___200630qox195223/output/C30A06597D_cfdna_N289-12_consensus.bait.pileup"
# pattern_bad <- ".*___[[:digit:]]{6}[a-z]{3}.*" #e.g. 210422140851sep
pattern <- ".*___[[:digit:]]{12}[a-z]{3}.*" #E.g. 200630qox195223
# 
# grepl(pattern, bad)
# grepl(pattern, good)

#These have the rigth pattern, i.e. the hg38 genome version (___)
hg38 <- cruk_files[ sapply(cruk_files, function(x) grepl(pattern = pattern, x = x) ) ]

#Grap the libraryid from these
hg38_libids <- sub(".*(N[[:digit:]]+-[[:digit:]]+).*", "\\1", hg38)


# Select a few samples from manifest that can be used for W3-4 development (These should be hg38)
#This is manifest with sample informations
manifest <- read.table("~/genomedk/matovanalysis/umiseq_analysis/R/cruk-manifest.csv", header = T, sep = " ", stringsAsFactors = F)

table(manifest$sample_type)
#adenom  control_no_cormorbiditet control_with_comorbiditet  control_with_IBD CRC_high_ctDNA CRC_pre-OP  other_cancer
#18                        18                        10               9               9           136       9 

#Add the hg version of the BAM files on the manifest file
manifest$hg38 <- sub(",.*$", "", manifest$library_id) %in% hg38_libids

#Which samples types have the hg38 / hg19 version 
table(manifest$hg38, manifest$sample_type)

#These are, say, high ctDNA samples (expected high signal) AND gh38 is TRUE - lets take 4 randomly
set.seed(8)
case <- sample(manifest[ manifest$sample_type == "CRC_high_ctDNA" & manifest$hg38  , "library_id"], 4) #"N289-224,N289-225" "N289-228" "N289-226,N289-227" "N289-221,N289-222"
control <- sample(manifest[ manifest$sample_type == "control_no_cormorbiditet" & manifest$hg38, "library_id"], 4) #"N289-115" "N289-174" "N289-269" "N289-282"


#Since some samples originates from two libraries ("split-samples") lest just take the first which is used for naming
case <- sub(",.*$", "", case)
control <- sub(",.*$", "", control)


#The frag pos files
fragpos_files <- list.files(pattern = "_consensus.txt") #the files fragpos
#Get the libraryids and convert e.g. N289_123 to N289-123  
fragpos_libids <- sub("_", "-", sub(".*(N[[:digit:]]+_[[:digit:]]+).*", "\\1", fragpos_files))


#Lets find those that we have selected as cases and controls
i <- sapply(c(case, control), function(x) grep(x, fragpos_libids) )

data <- 
  lapply(fragpos_files[i], function(x) {
  read.table(x, sep = "\t", header = F)[, -c(1:3)]
})
str(data, max.level = 1) #Its a list of pos-base * fraglen dataframes
#Make the list into an 3d array
data <- abind::abind(data, along = 3)
str(data, max.level = 1) #array of pos-base * fraglen * sample

dim1 <-  read.table(fragpos_files[1], sep = "\t", header = F)[, 1:3]
dim1 <- apply(dim1, 1, paste0, collapse = "_")
dim3 <- c( paste0(case, "_CRC"), paste0(control, "_CTL"))

dimnames(data) <- c(list(dim1, NULL, dim3))
saveRDS(data,  "test_fragpos.RDS") #Save the 4+4 sample test array of frag pos 
data <- readRDS("test_fragpos.RDS")  #load the test array

str(data)

#The first  case frag pos data
data[,,1] #A matrix - each row a substition

#The first  normal/ctrl frag pos data
data[,,5][1:4,150:175]


#Which are references  (those with many reads)
positions <- read.table("~/genomedk/PolyA/faststorage/BACKUP/CRUK/references/200504_sporacrc-coordinates-hg3819.tab", header = T)
head(positions)


raw <- read.table("C118A05812D_cfdna_N289_232_consensus.txt", sep = "\t", header = F)

raw[1:5,-c(1:3)]

getwd()


#raw <- read.table("sw620_control_N227-2640_consensus.txt", sep = "\t", header = F) #Read data
raw[1:10, 1:10] #Take a look at a corner of the data - first 3 columns are annotation (chr, pos, base)
hist(apply(raw[-c(1:3)], 1, sum), xlim = c(1500, 5000), ylim = c(0, 10000)) #Plot, excluding 0 counts.
