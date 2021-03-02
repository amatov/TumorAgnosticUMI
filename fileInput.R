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

cancer_SNPs

pt_id library_id sample_label cancer sample_type   pstage tumor_maf        sitemut_hg38

31   6329     N289-2   C22A06329D      1  CRC pre-OP pT2pN0M0     0.356  chr5:112815487_A/G

32   6329     N289-2   C22A06329D      1  CRC pre-OP pT2pN0M0     0.360  chr5:112838934_C/T

33   6329     N289-2   C22A06329D      1  CRC pre-OP pT2pN0M0     0.781   chr17:7673776_G/A

35   6330     N289-3   C47A06330D      1  CRC pre-OP pT1pN0M0     0.378  chr12:25245351_C/T

36   6399     N289-5   C27A06399D      1  CRC pre-OP pT2pN0M0     0.295  chr5:112792494_C/T

39   6414     N289-6   C42A06414D      1  CRC pre-OP pT2pN0M0     0.687  chr5:112792446_C/T

40   6414     N289-6   C42A06414D      1  CRC pre-OP pT2pN0M0     0.378  chr12:25245350_C/T

41   6485     N289-7   C32A06485D      1  CRC pre-OP pT2pN0M0     0.224  chr3:179218303_G/A

42   6485     N289-7   C32A06485D      1  CRC pre-OP pT2pN0M0     0.223  chr5:112815507_C/T

44   6485     N289-7   C32A06485D      1  CRC pre-OP pT2pN0M0     0.180   chr17:7674885_C/T

45   6485     N289-7   C32A06485D      1  CRC pre-OP pT2pN0M0     0.172   chr17:7674894_G/C

47   6543     N289-8   C26A06543D      1  CRC pre-OP pT2pN0M0     0.162 chr10:113165569_G/A

48   6543     N289-8   C26A06543D      1  CRC pre-OP pT2pN0M0     0.569  chr12:25245350_C/T

49   6543     N289-8   C26A06543D      1  CRC pre-OP pT2pN0M0     0.171  chr18:51078379_G/C

50   6551     N289-9   C40A06551D      1  CRC pre-OP pT2pN0M0     0.550  chr5:112815507_C/T

51   6551     N289-9   C40A06551D      1  CRC pre-OP pT2pN0M0     0.282  chr12:25245347_C/T

52   6551     N289-9   C40A06551D      1  CRC pre-OP pT2pN0M0     0.386   chr17:7673776_G/A

53   6551     N289-9   C40A06551D      1  CRC pre-OP pT2pN0M0     0.423   chr17:7675085_C/A

60   6610    N289-15   C28A06610D      1  CRC pre-OP pT2pN0M0     0.351  chr5:112819300_G/A

61   6610    N289-15   C28A06610D      1  CRC pre-OP pT2pN0M0     0.072  chr5:112839879_C/T

62   6610    N289-15   C28A06610D      1  CRC pre-OP pT2pN0M0     0.350  chr12:25245350_C/A

63   6610    N289-15   C28A06610D      1  CRC pre-OP pT2pN0M0     0.328   chr17:7674230_C/T

64   6610    N289-15   C28A06610D      1  CRC pre-OP pT2pN0M0     0.392  chr17:72123983_C/T

67   6625    N289-19   C31A06625D      1  CRC pre-OP pT2pN0Mx     0.686  chr5:112839576_C/T

68   6625    N289-19   C31A06625D      1  CRC pre-OP pT2pN0Mx     0.493  chr12:25225713_T/G

69   6625    N289-19   C31A06625D      1  CRC pre-OP pT2pN0Mx     0.670   chr17:7675161_G/A

70   6632    N289-20   C31A06632D      1  CRC pre-OP pT2pN0M0     0.180  chr4:152326215_G/A

71   6632    N289-20   C31A06632D      1  CRC pre-OP pT2pN0M0     0.240  chr5:112792446_C/T

73   6632    N289-20   C31A06632D      1  CRC pre-OP pT2pN0M0     0.264  chr12:25245347_C/T

74   6632    N289-20   C31A06632D      1  CRC pre-OP pT2pN0M0     0.548   chr17:7675088_C/T

75   6678    N289-21   C28A06678D      1  CRC pre-OP pT2pN0M0     0.556  chr5:112815487_A/G

76   6678    N289-21   C28A06678D      1  CRC pre-OP pT2pN0M0     0.355   chr17:7674256_T/C

78   6689    N289-22   C35A06689D      1  CRC pre-OP pT2pN0M0     0.119  chr4:152328233_G/A

79   6689    N289-22   C35A06689D      1  CRC pre-OP pT2pN0M0     0.261  chr5:112780895_C/T

80   6689    N289-22   C35A06689D      1  CRC pre-OP pT2pN0M0     0.179  chr5:112839501_C/T

81   6689    N289-22   C35A06689D      1  CRC pre-OP pT2pN0M0     0.149  chr12:25245328_C/A

82   6689    N289-22   C35A06689D      1  CRC pre-OP pT2pN0M0     0.169   chr17:7674953_T/C

83   6713    N289-23   C50A06713D      1  CRC pre-OP pT4pN0M0     0.146  chr5:112780895_C/T

85   6713    N289-23   C50A06713D      1  CRC pre-OP pT4pN0M0     0.143   chr17:7674221_G/A

86   6719    N289-24   C68A06719D      1  CRC pre-OP pT3pN2M0     0.083  chr1:118101793_T/C

87   6719    N289-24   C68A06719D      1  CRC pre-OP pT3pN2M0     0.671  chr5:112827248_G/C

88   6719    N289-24   C68A06719D      1  CRC pre-OP pT3pN2M0     0.333  chr5:112839450_G/T

89   6719    N289-24   C68A06719D      1  CRC pre-OP pT3pN2M0     0.692   chr17:7676044_A/T

90   6721    N289-25   C23A06721D      1  CRC pre-OP pT2pN0Mx     0.195  chr5:112815487_A/G

91   6721    N289-25   C23A06721D      1  CRC pre-OP pT2pN0Mx     0.268   chr17:7674917_T/C

92   6721    N289-25   C23A06721D      1  CRC pre-OP pT2pN0Mx     0.238  chr17:72123748_C/A

93   6726    N289-26   C31A06726D      1  CRC pre-OP pT2pN1Mx     0.650  chr5:112839729_G/T

94   6726    N289-26   C31A06726D      1  CRC pre-OP pT2pN1Mx     0.283  chr12:25245347_C/T

98   6759    N289-29   C29A06759D      1  CRC pre-OP pT2pN0M0     0.308  chr5:112838741_G/A

100  6759    N289-29   C29A06759D      1  CRC pre-OP pT2pN0M0     0.667  chr12:25227341_T/A

101  6759    N289-29   C29A06759D      1  CRC pre-OP pT2pN0M0     0.403   chr17:7675088_C/T

102  6766    N289-30   C66A06766D      1  CRC pre-OP pT4pN0M0     0.176  chr3:179218303_G/A

103  6766    N289-30   C66A06766D      1  CRC pre-OP pT4pN0M0     0.347  chr5:112780895_C/T

104  6766    N289-30   C66A06766D      1  CRC pre-OP pT4pN0M0     0.495  chr5:112815487_A/G

105  6766    N289-30   C66A06766D      1  CRC pre-OP pT4pN0M0     0.451  chr12:25245350_C/T

106  6780    N289-31   C36A06780D      1  CRC pre-OP pT2pN0M0     0.484  chr5:112838196_G/T

107  6780    N289-31   C36A06780D      1  CRC pre-OP pT2pN0M0     0.121  chr7:140753355_C/T

108  6780    N289-31   C36A06780D      1  CRC pre-OP pT2pN0M0     0.254   chr17:7673796_C/T

111  6792    N289-33   C62A06792D      1  CRC pre-OP pT3pN1M0     0.237  chr5:112815487_A/G

112  6792    N289-33   C62A06792D      1  CRC pre-OP pT3pN1M0     0.348   chr17:7674250_C/T

114  6815    N289-35   C67A06815D      1  CRC pre-OP pT3pN2M0     0.218  chr5:112838116_T/G

115  6815    N289-35   C67A06815D      1  CRC pre-OP pT3pN2M0     0.254  chr12:25245347_C/T

116  6815    N289-35   C67A06815D      1  CRC pre-OP pT3pN2M0     0.436   chr17:7674221_G/A

117  6815    N289-35   C67A06815D      1  CRC pre-OP pT3pN2M0     0.204   chrX:64191396_G/A

118  6823    N289-36   C28A06823D      1  CRC pre-OP pT1pN0M0     0.519  chr5:112838740_G/A

119  6823    N289-36   C28A06823D      1  CRC pre-OP pT1pN0M0     0.303   chr17:7674220_C/T

120  6825    N289-37   C42A06825D      1  CRC pre-OP pT2pN0M0     0.369  chr5:112839427_C/A

122  6825    N289-37   C42A06825D      1  CRC pre-OP pT2pN0M0     0.405   chr17:7675993_C/A

125  6901    N289-41   C54A06901D      1  CRC pre-OP pT3pN2M0     0.290   chr17:7675127_A/C

129  6917    N289-45   C54A06917D      1  CRC pre-OP pT3pN2M0     0.337  chr12:25245347_C/T

130  6917    N289-45   C54A06917D      1  CRC pre-OP pT3pN2M0     0.582   chr17:7675088_C/T

131  6917    N289-45   C54A06917D      1  CRC pre-OP pT3pN2M0     0.029   chrX:64191300_C/G

132  6918    N289-46   C66A06918D      1  CRC pre-OP pT3pN2M0     0.143  chr7:140753336_A/T

133  6918    N289-46   C66A06918D      1  CRC pre-OP pT3pN2M0     0.154   chr17:7674953_T/C

134  6918    N289-46   C66A06918D      1  CRC pre-OP pT3pN2M0     0.138   chrX:64191455_C/A

137  6942    N289-84   C20A06942D      1  CRC pre-OP pT2pN0M0     0.307  chr12:25245350_C/A

138  6948    N289-83   C58A06948D      1  CRC pre-OP pT3pN1M0     0.091  chr3:179218304_A/G

139  6948    N289-83   C58A06948D      1  CRC pre-OP pT3pN1M0     0.122  chr12:25245350_C/T

141  6966    N289-50   C67A06966D      1  CRC pre-OP pT2pN1M0     0.685  chr5:112839810_C/T

142  6966    N289-50   C67A06966D      1  CRC pre-OP pT2pN1M0     0.655  chr12:25227342_T/A

144  6969    N289-85   C44A06969D      1  CRC pre-OP pT2pN0M0     0.382  chr5:112792494_C/T

145  6969    N289-85   C44A06969D      1  CRC pre-OP pT2pN0M0     0.544   chr17:7675088_C/T

146  6975    N289-52   C41A06975D      1  CRC pre-OP pT2pN1M0     0.433  chr5:112839661_C/G

147  6975    N289-52   C41A06975D      1  CRC pre-OP pT2pN1M0     0.345   chr17:7670685_G/A

148  6975    N289-52   C41A06975D      1  CRC pre-OP pT2pN1M0     0.322   chr17:7673776_G/A

149  6996    N289-53   C65A06996D      1  CRC pre-OP pT3pN1M0     0.171  chr1:114713908_T/C

150  6996    N289-53   C65A06996D      1  CRC pre-OP pT3pN1M0     0.166  chr5:112819266_C/T

153  7005    N289-56   C17A07005D      1  CRC pre-OP pT2pN0M0     0.461  chr4:152328233_G/A

154  7005    N289-56   C17A07005D      1  CRC pre-OP pT2pN0M0     0.614  chr12:25245347_C/T

155  7005    N289-56   C17A07005D      1  CRC pre-OP pT2pN0M0     0.810   chr17:7673609_C/T

156  7005    N289-56   C17A07005D      1  CRC pre-OP pT2pN0M0     0.825  chr18:51078378_T/A

157  7006    N289-57   C53A07006D      1  CRC pre-OP pT3pN1M0     0.263  chr3:179199066_G/A

158  7006    N289-57   C53A07006D      1  CRC pre-OP pT3pN1M0     0.467  chr12:25245350_C/A

159  7007    N289-58   C47A07007D      1  CRC pre-OP pT3pN1M0     0.269  chr5:112838934_C/T

161  7007    N289-58   C47A07007D      1  CRC pre-OP pT3pN1M0     0.113   chr17:7674221_G/C

163  7029    N289-61   C55A07029D      1  CRC pre-OP pT2pN0M0     0.375  chr1:114713908_T/A

164  7029    N289-61   C55A07029D      1  CRC pre-OP pT2pN0M0     0.672  chr5:112839510_G/T

165  7029    N289-61   C55A07029D      1  CRC pre-OP pT2pN0M0     0.337   chr17:7673776_G/A

168  7037    N289-63   C43A07037D      1  CRC pre-OP pT2pN0M0     0.439  chr5:112838455_T/A

169  7037    N289-63   C43A07037D      1  CRC pre-OP pT2pN0M0     0.250  chr5:112839879_C/T