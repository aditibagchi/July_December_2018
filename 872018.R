primary_medulloblastoma_samples_maf <- read.delim("~/Desktop/Data_Analysis/Data_Files/Exome_Data/primary_medulloblastoma_samples_maf.txt", header=FALSE, comment.char="#")
maf_1 <- read.maf(maf = "~/Desktop/Data_Analysis/Data_Files/Exome_Data/primary_medulloblastoma_samples_maf.txt" )
mafSummary(maf_1)
subsetMaf(maf_1, tsb = "ABSJ052018CHTN-7-PM")

nrow(ABSJ052018CHTN7PM)

ncol(ABSJ052018CHTN7PM)
ABSJ052018CHTN7PM_DF <- as.data.frame(ABSJ052018CHTN7PM)
keep <- 3:(ncol(ABSJ052018CHTN7PM_DF) - 113)

ABSJ052018CHTN7PM_1 <- ABSJ052018CHTN7PM_DF[, keep]
head(ABSJ052018CHTN7PM_1)

maf_2 <- read.maf(maf = "~/Desktop/Data_Analysis/Data_Files/Exome_Data/recurrent_medulloblastoma_samples_maf.txt" )
mafSummary(maf_2)
subsetMaf(maf_1, tsb = "ABSJ052018CHTN-7-RM")
subsetMaf(maf_2, tsb = "ABSJ052018CHTN-7-RM")
ABSJ052018CHTN7RM <- subsetMaf(maf_2, tsb = "ABSJ052018CHTN-7-RM")
nrow(ABSJ052018CHTN7RM)
NCOL(ABSJ052018CHTN7RM)
ABSJ052018CHTN7RM_DF <- as.data.frame(ABSJ052018CHTN7RM)
keep1 <- 3:ncol(ABSJ052018CHTN7RM_DF) - 113)
keep1 <- 3:(ncol(ABSJ052018CHTN7RM_DF) - 113)
ABSJ052018CHTN7RM_1 <- ABSJ052018CHTN7RM_DF[, keep1]

