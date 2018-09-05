## Read the maf file using maftools               
maf1 <- read.maf(maf = "~/Desktop/Data_Analysis/Data_Files/Exome_Data/primary_medulloblastoma_samples_maf.txt"
mafSummary(maf1)
plotmafSummary(maf1)


## creating a data frame and then seting a minimumdepth filter of 100X 

Maf1_MB11 <- subsetMaf(maf1, tsb = "ABSJ112017MB11"  , genes = NULL, fields = NULL, query = NULL,
                      mafObj = FALSE, includeSyn = TRUE, isTCGA = FALSE, restrictTo = "all")
str(Maf1_MB11)
min.depth.filter  = 100
Maf1_MB11_depth = subset(Maf1_MB11,Maf1_MB11$t_depth>=min.depth.filter)
str(Maf1_MB11_depth)
cat('Maf1_MB11_depth in dbSNP: ',length(which(! Maf1_MB11_depth$dbSNP_RS=='novel')))
cat('Maf1_MB11_depth NOT dbSNP: ',length(which(Maf1_MB11_depth$dbSNP_RS=='novel')))

filter = 3
Maf1_MB11_depth_AF3 = subset(Maf1_MB11_depth,Maf1_MB11_depth$AF>=filter)
str(Maf1_MB11_depth_AF3)

filter = 1
Maf1_MB11_depth_AF1 = subset(Maf1_MB11_depth,Maf1_MB11_depth$AF>=filter)
str(Maf1_MB11_depth_AF1)

filter = 0.5
Maf1_MB11_depth_AF = subset(Maf1_MB11_depth,Maf1_MB11_depth$AF>=filter)
str(Maf1_MB11_depth_AF)
cat('Maf1_MB11_depth_AF in dbSNP: ',length(which(! Maf1_MB11_depth_AF$dbSNP_RS=='novel')))
cat('Maf1_MB11_depth_AF NOT dbSNP: ',length(which(Maf1_MB11_depth_AF$dbSNP_RS=='novel')))




View(Maf1_df)
Maf1_BR2 <- subsetMaf(maf1, tsb = "ABSJ112017BR2Comb"  , genes = NULL, fields = NULL, query = NULL,
                     mafObj = FALSE, includeSyn = TRUE, isTCGA = FALSE, restrictTo = "all")
str(Maf1_BR2)
min.depth.filter  = 100
Maf1_BR2_depth = subset(Maf1_BR2,Maf1_BR2$t_depth>=min.depth.filter)
str(Maf1_BR2_depth)
cat('Maf1_BR2_depth in dbSNP: ',length(which(! Maf1_BR2_depth$dbSNP_RS=='novel')))
cat('Maf1_BR2_depth NOT dbSNP: ',length(which(Maf1_BR2_depth$dbSNP_RS=='novel')))





Maf1_MB2 <- subsetMaf(maf1, tsb = "ABSJ112017MB2"  , genes = NULL, fields = NULL, query = NULL,
                      mafObj = FALSE, includeSyn = TRUE, isTCGA = FALSE, restrictTo = "all")
str(Maf1_MB2)
min.depth.filter  = 100
Maf1_MB2_depth = subset(Maf1_MB2,Maf1_MB2$t_depth>=min.depth.filter)
str(Maf1_MB2_depth)
cat('Maf1_MB2_depth in dbSNP: ',length(which(! Maf1_MB2_depth$dbSNP_RS=='novel')))
cat('Maf1_MB2_depth NOT dbSNP: ',length(which(Maf1_MB2_depth$dbSNP_RS=='novel')))
filter = 3
Maf1_MB2_depth_AF3 = subset(Maf1_MB2_depth,Maf1_MB2_depth$AF>=filter)
str(Maf1_MB2_depth_AF3)

filter = 1
Maf1_MB2_depth_AF1 = subset(Maf1_MB2_depth,Maf1_MB2_depth$AF>=filter)
str(Maf1_MB2_depth_AF1)

filter = 0.5
Maf1_MB2_depth_AF = subset(Maf1_MB2_depth,Maf1_MB2_depth$AF>=filter)
str(Maf1_MB2_depth_AF)
cat('Maf1_MB2_depth_AF in dbSNP: ',length(which(! Maf1_MB2_depth_AF$dbSNP_RS=='novel')))
cat('Maf1_MB2_depth_AF NOT dbSNP: ',length(which(Maf1_MB2_depth_AF$dbSNP_RS=='novel')))





Maf1_MB28 <- subsetMaf(maf1, tsb = "ABSJ112017MB28"  , genes = NULL, fields = NULL, query = NULL,
                      mafObj = FALSE, includeSyn = TRUE, isTCGA = FALSE, restrictTo = "all")
str(Maf1_MB28)
min.depth.filter  = 100
Maf1_MB28_depth = subset(Maf1_MB28,Maf1_MB28$t_depth>=min.depth.filter)
str(Maf1_MB28_depth)

cat('Maf1_MB28_depth in dbSNP: ',length(which(! Maf1_MB28_depth$dbSNP_RS=='novel')))
cat('Maf1_MB28_depth NOT dbSNP: ',length(which(Maf1_MB28_depth$dbSNP_RS=='novel')))

filter = 3
Maf1_MB28_depth_AF3 = subset(Maf1_MB28_depth,Maf1_MB28_depth$AF>=filter)
str(Maf1_MB28_depth_AF3)

filter = 1
Maf1_MB28_depth_AF1 = subset(Maf1_MB28_depth,Maf1_MB28_depth$AF>=filter)
str(Maf1_MB28_depth_AF1)

filter = 0.5
Maf1_MB28_depth_AF = subset(Maf1_MB28_depth,Maf1_MB28_depth$AF>=filter)
str(Maf1_MB28_depth_AF)
cat('Maf1_MB28_depth_AF in dbSNP: ',length(which(! Maf1_MB28_depth_AF$dbSNP_RS=='novel')))
cat('Maf1_MB28_depth_AF NOT dbSNP: ',length(which(Maf1_MB28_depth_AF$dbSNP_RS=='novel')))


Maf1_MB13 <- subsetMaf(maf1, tsb = "ABSJ112017MB13"  , genes = NULL, fields = NULL, query = NULL,
                       mafObj = FALSE, includeSyn = TRUE, isTCGA = FALSE, restrictTo = "all")
min.depth.filter  = 100
Maf1_MB13_depth = subset(Maf1_MB13,Maf1_MB13$t_depth>=min.depth.filter)
str(Maf1_MB13_depth)

cat('Maf1_MB13_depth in dbSNP: ',length(which(! Maf1_MB13_depth$dbSNP_RS=='novel')))
cat('Maf1_MB13_depth NOT dbSNP: ',length(which(Maf1_MB13_depth$dbSNP_RS=='novel')))

filter = 3
Maf1_MB13_depth_AF3 = subset(Maf1_MB13_depth,Maf1_MB13_depth$AF>=filter)
str(Maf1_MB13_depth_AF3)

filter = 1
Maf1_MB13_depth_AF1 = subset(Maf1_MB13_depth,Maf1_MB13_depth$AF>=filter)
str(Maf1_MB13_depth_AF1)

filter = 0.5
Maf1_MB13_depth_AF = subset(Maf1_MB13_depth,Maf1_MB13_depth$AF>=filter)
str(Maf1_MB13_depth_AF)
cat('Maf1_MB13_depth_AF in dbSNP: ',length(which(! Maf1_MB13_depth_AF$dbSNP_RS=='novel')))
cat('Maf1_MB13_depth_AF NOT dbSNP: ',length(which(Maf1_MB13_depth_AF$dbSNP_RS=='novel')))


Maf1_MB9 <- subsetMaf(maf1, tsb = "ABSJ112017MB9"  , genes = NULL, fields = NULL, query = NULL,
                       mafObj = FALSE, includeSyn = TRUE, isTCGA = FALSE, restrictTo = "all")
min.depth.filter  = 100
Maf1_MB9_depth = subset(Maf1_MB9,Maf1_MB9$t_depth>=min.depth.filter)
str(Maf1_MB9_depth)

cat('Maf1_MB9_depth in dbSNP: ',length(which(! Maf1_MB9_depth$dbSNP_RS=='novel')))
cat('Maf1_MB9_depth NOT dbSNP: ',length(which(Maf1_MB9_depth$dbSNP_RS=='novel')))

filter = 3
Maf1_MB9_depth_AF3 = subset(Maf1_MB9_depth,Maf1_MB9_depth$AF>=filter)
str(Maf1_MB13_depth_AF3)

filter = 1
Maf1_MB9_depth_AF1 = subset(Maf1_MB9_depth,Maf1_MB9_depth$AF>=filter)
str(Maf1_MB9_depth_AF1)

filter = 0.5
Maf1_MB9_depth_AF = subset(Maf1_MB9_depth,Maf1_MB9_depth$AF>=filter)
str(Maf1_MB9_depth_AF)
cat('Maf1_MB9_depth_AF in dbSNP: ',length(which(! Maf1_MB9_depth_AF$dbSNP_RS=='novel')))
cat('Maf1_MB9_depth_AF NOT dbSNP: ',length(which(Maf1_MB9_depth_AF$dbSNP_RS=='novel')))

Maf1_MB7 <- subsetMaf(maf1, tsb = "ABSJ112017MB7"  , genes = NULL, fields = NULL, query = NULL,
                      mafObj = FALSE, includeSyn = TRUE, isTCGA = FALSE, restrictTo = "all")
min.depth.filter  = 100
Maf1_MB7_depth = subset(Maf1_MB7,Maf1_MB7$t_depth>=min.depth.filter)
str(Maf1_MB7_depth)

cat('Maf1_MB7_depth in dbSNP: ',length(which(! Maf1_MB7_depth$dbSNP_RS=='novel')))
cat('Maf1_MB7_depth NOT dbSNP: ',length(which(Maf1_MB7_depth$dbSNP_RS=='novel')))


filter = 3
Maf1_MB7_depth_AF3 = subset(Maf1_MB7_depth,Maf1_MB7_depth$AF>=filter)
str(Maf1_MB7_depth_AF3)

filter = 1
Maf1_MB7_depth_AF1 = subset(Maf1_MB7_depth,Maf1_MB7_depth$AF>=filter)
str(Maf1_MB7_depth_AF1)

filter = 0.5
Maf1_MB7_depth_AF = subset(Maf1_MB7_depth,Maf1_MB7_depth$AF>=filter)
str(Maf1_MB7_depth_AF)
cat('Maf1_MB7_depth_AF in dbSNP: ',length(which(! Maf1_MB7_depth_AF$dbSNP_RS=='novel')))
cat('Maf1_MB7_depth_AF NOT dbSNP: ',length(which(Maf1_MB7_depth_AF$dbSNP_RS=='novel')))



Maf1_MB14 <- subsetMaf(maf1, tsb = "ABSJ112017MB14"  , genes = NULL, fields = NULL, query = NULL,
                     mafObj = FALSE, includeSyn = TRUE, isTCGA = FALSE, restrictTo = "all")
min.depth.filter  = 100
Maf1_MB14_depth = subset(Maf1_MB14,Maf1_MB14$t_depth>=min.depth.filter)
str(Maf1_MB14_depth)

cat('Maf1_MB14_depth in dbSNP: ',length(which(! Maf1_MB14_depth$dbSNP_RS=='novel')))
cat('Maf1_MB14_depth NOT dbSNP: ',length(which(Maf1_MB14_depth$dbSNP_RS=='novel')))

filter = 3
Maf1_MB14_depth_AF3 = subset(Maf1_MB14_depth,Maf1_MB14_depth$AF>=filter)
str(Maf1_MB14_depth_AF3)

filter = 1
Maf1_MB14_depth_AF1 = subset(Maf1_MB14_depth,Maf1_MB14_depth$AF>=filter)
str(Maf1_MB14_depth_AF1)

filter = 0.5
Maf1_MB14_depth_AF = subset(Maf1_MB14_depth,Maf1_MB14_depth$AF>=filter)
str(Maf1_MB14_depth_AF)
cat('Maf1_MB14_depth_AF in dbSNP: ',length(which(! Maf1_MB14_depth_AF$dbSNP_RS=='novel')))
cat('Maf1_MB14_depth_AF NOT dbSNP: ',length(which(Maf1_MB14_depth_AF$dbSNP_RS=='novel')))



Maf1_MB27 <- subsetMaf(maf1, tsb = "ABSJ112017MB27"  , genes = NULL, fields = NULL, query = NULL,
                       mafObj = FALSE, includeSyn = TRUE, isTCGA = FALSE, restrictTo = "all")
min.depth.filter  = 100
Maf1_MB27_depth = subset(Maf1_MB27,Maf1_MB27$t_depth>=min.depth.filter)
str(Maf1_MB27_depth)

cat('Maf1_MB27_depth in dbSNP: ',length(which(! Maf1_MB27_depth$dbSNP_RS=='novel')))
cat('Maf1_MB27_depth NOT dbSNP: ',length(which(Maf1_MB27_depth$dbSNP_RS=='novel')))

filter = 3
Maf1_MB27_depth_AF3 = subset(Maf1_MB27_depth,Maf1_MB27_depth$AF>=filter)
str(Maf1_MB27_depth_AF3)

filter = 1
Maf1_MB27_depth_AF1 = subset(Maf1_MB27_depth,Maf1_MB27_depth$AF>=filter)
str(Maf1_MB27_depth_AF1)

filter = 0.5
Maf1_MB27_depth_AF = subset(Maf1_MB27_depth,Maf1_MB27_depth$AF>=filter)
str(Maf1_MB27_depth_AF)
cat('Maf1_MB27_depth_AF in dbSNP: ',length(which(! Maf1_MB27_depth_AF$dbSNP_RS=='novel')))
cat('Maf1_MB27_depth_AF NOT dbSNP: ',length(which(Maf1_MB27_depth_AF$dbSNP_RS=='novel')))





Maf1_MB22 <- subsetMaf(maf1, tsb = "ABSJ112017MB22"  , genes = NULL, fields = NULL, query = NULL,
                       mafObj = FALSE, includeSyn = TRUE, isTCGA = FALSE, restrictTo = "all")
min.depth.filter  = 100
Maf1_MB22_depth = subset(Maf1_MB22,Maf1_MB22$t_depth>=min.depth.filter)
str(Maf1_MB22_depth)

cat('Maf1_MB22_depth in dbSNP: ',length(which(! Maf1_MB22_depth$dbSNP_RS=='novel')))
cat('Maf1_MB22_depth NOT dbSNP: ',length(which(Maf1_MB22_depth$dbSNP_RS=='novel')))

filter = 3
Maf1_MB22_depth_AF3 = subset(Maf1_MB22_depth,Maf1_MB22_depth$AF>=filter)
str(Maf1_MB22_depth_AF3)

filter = 1
Maf1_MB22_depth_AF1 = subset(Maf1_MB22_depth,Maf1_MB22_depth$AF>=filter)
str(Maf1_MB22_depth_AF1)

filter = 0.5
Maf1_MB22_depth_AF = subset(Maf1_MB22_depth,Maf1_MB22_depth$AF>=filter)
str(Maf1_MB22_depth_AF)
cat('Maf1_MB22_depth_AF in dbSNP: ',length(which(! Maf1_MB22_depth_AF$dbSNP_RS=='novel')))
cat('Maf1_MB22_depth_AF NOT dbSNP: ',length(which(Maf1_MB22_depth_AF$dbSNP_RS=='novel')))



Maf1_MB5 <- subsetMaf(maf1, tsb = "ABSJ112017MB5"  , genes = NULL, fields = NULL, query = NULL,
                       mafObj = FALSE, includeSyn = TRUE, isTCGA = FALSE, restrictTo = "all")
min.depth.filter  = 100
Maf1_MB5_depth = subset(Maf1_MB5,Maf1_MB5$t_depth>=min.depth.filter)
str(Maf1_MB5_depth)

cat('Maf1_MB5_depth in dbSNP: ',length(which(! Maf1_MB5_depth$dbSNP_RS=='novel')))
cat('Maf1_MB5_depth NOT dbSNP: ',length(which(Maf1_MB5_depth$dbSNP_RS=='novel')))

filter = 3
Maf1_MB5_depth_AF3 = subset(Maf1_MB5_depth,Maf1_MB5_depth$AF>=filter)
str(Maf1_MB5_depth_AF3)

filter = 1
Maf1_MB5_depth_AF1 = subset(Maf1_MB5_depth,Maf1_MB5_depth$AF>=filter)
str(Maf1_MB5_depth_AF1)

filter = 0.5
Maf1_MB5_depth_AF = subset(Maf1_MB5_depth,Maf1_MB5_depth$AF>=filter)
str(Maf1_MB5_depth_AF)
cat('Maf1_MB5_depth_AF in dbSNP: ',length(which(! Maf1_MB5_depth_AF$dbSNP_RS=='novel')))
cat('Maf1_MB5_depth_AF NOT dbSNP: ',length(which(Maf1_MB5_depth_AF$dbSNP_RS=='novel')))

Maf1_MB19 <- subsetMaf(maf1, tsb = "ABSJ112017MB19"  , genes = NULL, fields = NULL, query = NULL,
                      mafObj = FALSE, includeSyn = TRUE, isTCGA = FALSE, restrictTo = "all")
min.depth.filter  = 100
Maf1_MB19_depth = subset(Maf1_MB19,Maf1_MB19$t_depth>=min.depth.filter)
str(Maf1_MB19_depth)

cat('Maf1_MB19_depth in dbSNP: ',length(which(! Maf1_MB19_depth$dbSNP_RS=='novel')))
cat('Maf1_MB19_depth NOT dbSNP: ',length(which(Maf1_MB19_depth$dbSNP_RS=='novel')))


filter = 3
Maf1_MB19_depth_AF3 = subset(Maf1_MB19_depth,Maf1_MB19_depth$AF>=filter)
str(Maf1_MB19_depth_AF3)

filter = 1
Maf1_MB19_depth_AF1 = subset(Maf1_MB19_depth,Maf1_MB19_depth$AF>=filter)
str(Maf1_MB19_depth_AF1)

filter = 0.5
Maf1_MB19_depth_AF = subset(Maf1_MB19_depth,Maf1_MB19_depth$AF>=filter)
str(Maf1_MB19_depth_AF)
cat('Maf1_MB19_depth_AF in dbSNP: ',length(which(! Maf1_MB19_depth_AF$dbSNP_RS=='novel')))
cat('Maf1_MB19_depth_AF NOT dbSNP: ',length(which(Maf1_MB19_depth_AF$dbSNP_RS=='novel')))


Maf1_MB17 <- subsetMaf(maf1, tsb = "ABSJ112017MB17"  , genes = NULL, fields = NULL, query = NULL,
                       mafObj = FALSE, includeSyn = TRUE, isTCGA = FALSE, restrictTo = "all")
min.depth.filter  = 100
Maf1_MB17_depth = subset(Maf1_MB17,Maf1_MB17$t_depth>=min.depth.filter)
str(Maf1_MB17_depth)

cat('Maf1_MB17_depth in dbSNP: ',length(which(! Maf1_MB17_depth$dbSNP_RS=='novel')))
cat('Maf1_MB17_depth NOT dbSNP: ',length(which(Maf1_MB17_depth$dbSNP_RS=='novel')))

filter = 3
Maf1_MB17_depth_AF3 = subset(Maf1_MB17_depth,Maf1_MB17_depth$AF>=filter)
str(Maf1_MB17_depth_AF3)

filter = 1
Maf1_MB17_depth_AF1 = subset(Maf1_MB17_depth,Maf1_MB17_depth$AF>=filter)
str(Maf1_MB17_depth_AF1)

filter = 0.5
Maf1_MB17_depth_AF = subset(Maf1_MB17_depth,Maf1_MB5_depth$AF>=filter)
str(Maf1_MB17_depth_AF)
cat('Maf1_MB17_depth_AF in dbSNP: ',length(which(! Maf1_MB5_depth_AF$dbSNP_RS=='novel')))
cat('Maf1_MB17_depth_AF NOT dbSNP: ',length(which(Maf1_MB5_depth_AF$dbSNP_RS=='novel')))



Maf1_MB21 <- subsetMaf(maf1, tsb = "ABSJ112017MB21"  , genes = NULL, fields = NULL, query = NULL,
                       mafObj = FALSE, includeSyn = TRUE, isTCGA = FALSE, restrictTo = "all")
min.depth.filter  = 100
Maf1_MB21_depth = subset(Maf1_MB21,Maf1_MB21$t_depth>=min.depth.filter)
str(Maf1_MB21_depth)
cat('Maf1_MB21_depth in dbSNP: ',length(which(! Maf1_MB21_depth$dbSNP_RS=='novel')))
cat('Maf1_MB21_depth NOT dbSNP: ',length(which(Maf1_MB21_depth$dbSNP_RS=='novel'))) 

filter = 3
Maf1_MB21_depth_AF3 = subset(Maf1_MB21_depth,Maf1_MB21_depth$AF>=filter)
str(Maf1_MB21_depth_AF3)

filter = 1
Maf1_MB21_depth_AF1 = subset(Maf1_MB21_depth,Maf1_MB17_depth$AF>=filter)
str(Maf1_MB17_depth_AF1)

filter = 0.5
Maf1_MB21_depth_AF = subset(Maf1_MB21_depth,Maf1_MB21_depth$AF>=filter)
str(Maf1_MB21_depth_AF)
cat('Maf1_MB21_depth_AF in dbSNP: ',length(which(! Maf1_MB21_depth_AF$dbSNP_RS=='novel')))
cat('Maf1_MB21_depth_AF NOT dbSNP: ',length(which(Maf1_MB21_depth_AF$dbSNP_RS=='novel')))

Maf1_MB25 <- subsetMaf(maf1, tsb = "ABSJ112017MB25"  , genes = NULL, fields = NULL, query = NULL,
                       mafObj = FALSE, includeSyn = TRUE, isTCGA = FALSE, restrictTo = "all")
min.depth.filter  = 100
Maf1_MB25_depth = subset(Maf1_MB25,Maf1_MB25$t_depth>=min.depth.filter)
str(Maf1_MB25_depth)

cat('Maf1_MB25_depth in dbSNP: ',length(which(! Maf1_MB25_depth$dbSNP_RS=='novel')))
cat('Maf1_MB25_depth NOT dbSNP: ',length(which(Maf1_MB25_depth$dbSNP_RS=='novel'))) 

filter = 3
Maf1_MB25_depth_AF3 = subset(Maf1_MB25_depth,Maf1_MB25_depth$AF>=filter)
str(Maf1_MB25_depth_AF3)

filter = 1
Maf1_MB25_depth_AF1 = subset(Maf1_MB25_depth,Maf1_MB25_depth$AF>=filter)
str(Maf1_MB25_depth_AF1)

filter = 0.5
Maf1_MB25_depth_AF = subset(Maf1_MB25_depth,Maf1_MB25_depth$AF>=filter)
str(Maf1_MB25_depth_AF)
cat('Maf1_MB25_depth_AF in dbSNP: ',length(which(! Maf1_MB25_depth_AF$dbSNP_RS=='novel')))
cat('Maf1_MB25_depth_AF NOT dbSNP: ',length(which(Maf1_MB25_depth_AF$dbSNP_RS=='novel')))


Maf1_MB_CHTN_PM_3 <- subsetMaf(maf1, tsb = "ABSJ052018CHTN-3-PM"  , genes = NULL, fields = NULL, query = NULL,
                               +                        mafObj = FALSE, includeSyn = TRUE, isTCGA = FALSE, restrictTo = "all")

min.depth.filter  = 100
Maf1_MB_CHTN_PM_3_depth = subset(Maf1_MB_CHTN_PM_3,Maf1_MB_CHTN_PM_3$t_depth>=min.depth.filter)
str(Maf1_MB_CHTN_PM_3_depth)

cat('Maf1_MB_CHTN_PM_3_depth in dbSNP: ',length(which(! Maf1_MB_CHTN_PM_3_depth$dbSNP_RS=='novel')))
cat('Maf1_MB_CHTN_PM_3_depth NOT dbSNP: ',length(which(Maf1_MB_CHTN_PM_3_depth$dbSNP_RS=='novel'))) 



filter = 3
Maf1_MB_CHTN_PM_3_depth_AF3 = subset(Maf1_MB_CHTN_PM_3_depth,Maf1_MB_CHTN_PM_3_depth$AF>=filter)
str(Maf1_MB_CHTN_PM_3_depth_AF3)

filter = 1
Maf1_MB_CHTN_PM_3_depth_AF1 = subset(Maf1_MB_CHTN_PM_3_depth,Maf1_MB_CHTN_PM_3_depth$AF>=filter)
str(Maf1_MB_CHTN_PM_3_depth_AF1)

filter = 0.5
Maf1_MB_CHTN_PM_3_depth_AF = subset(Maf1_MB_CHTN_PM_3_depth,Maf1_MB_CHTN_PM_3_depth$AF>=filter)
str(Maf1_MB_CHTN_PM_3_depth_AF)
cat('Maf1_MB_CHTN_PM_3_depth_AF in dbSNP: ',length(which(! Maf1_MB_CHTN_PM_3_depth_AF$dbSNP_RS=='novel')))
cat(' Maf1_MB_CHTN_PM_3_depth_AF NOT dbSNP: ',length(which(Maf1_MB_CHTN_PM_3_depth_AF$dbSNP_RS=='novel')))

Maf1_MB_CHTN_PM_5 <- subsetMaf(maf1, tsb = "ABSJ052018CHTN-5-PM"  , genes = NULL, fields = NULL, query = NULL,
                               mafObj = FALSE, includeSyn = TRUE, isTCGA = FALSE, restrictTo = "all")

min.depth.filter  = 100
Maf1_MB_CHTN_PM_5_depth = subset(Maf1_MB_CHTN_PM_5,Maf1_MB_CHTN_PM_5$t_depth>=min.depth.filter)
str(Maf1_MB_CHTN_PM_5_depth)

cat('Maf1_MB_CHTN_PM_5_depth in dbSNP: ',length(which(! Maf1_MB_CHTN_PM_5_depth$dbSNP_RS=='novel')))
cat('Maf1_MB_CHTN_PM_5_depth NOT dbSNP: ',length(which(Maf1_MB_CHTN_PM_5_depth$dbSNP_RS=='novel')))

filter = 3
Maf1_MB_CHTN_PM_5_depth_AF3 = subset(Maf1_MB_CHTN_PM_5_depth,Maf1_MB_CHTN_PM_5_depth$AF>=filter)
str(Maf1_MB_CHTN_PM_5_depth_AF3)

filter = 1
Maf1_MB_CHTN_PM_5_depth_AF1 = subset(Maf1_MB_CHTN_PM_5_depth,Maf1_MB_CHTN_PM_5_depth$AF>=filter)
str(Maf1_MB_CHTN_PM_5_depth_AF1)

filter = 0.5
Maf1_MB_CHTN_PM_5_depth_AF = subset(Maf1_MB_CHTN_PM_5_depth,Maf1_MB_CHTN_PM_5_depth$AF>=filter)
str(Maf1_MB_CHTN_PM_5_depth_AF)
cat('Maf1_MB_CHTN_PM_5_depth_AF in dbSNP: ',length(which(! Maf1_MB_CHTN_PM_5_depth_AF$dbSNP_RS=='novel')))
cat(' Maf1_MB_CHTN_PM_5_depth_AF NOT dbSNP: ',length(which(Maf1_MB_CHTN_PM_5_depth_AF$dbSNP_RS=='novel')))

Maf1_MB_CHTN_PM_2 <- subsetMaf(maf1, tsb = "ABSJ052018CHTN-2-PM"  , genes = NULL, fields = NULL, query = NULL,
                               mafObj = FALSE, includeSyn = TRUE, isTCGA = FALSE, restrictTo = "all")

min.depth.filter  = 100
Maf1_MB_CHTN_PM_2_depth = subset(Maf1_MB_CHTN_PM_2,Maf1_MB_CHTN_PM_2$t_depth>=min.depth.filter)
str(Maf1_MB_CHTN_PM_2_depth)
cat('Maf1_MB_CHTN_PM_2_depth in dbSNP: ',length(which(! Maf1_MB_CHTN_PM_2_depth$dbSNP_RS=='novel')))
cat('Maf1_MB_CHTN_PM_2_depth NOT dbSNP: ',length(which(Maf1_MB_CHTN_PM_2_depth$dbSNP_RS=='novel')))

filter = 3
Maf1_MB_CHTN_PM_2_depth_AF3 = subset(Maf1_MB_CHTN_PM_2_depth,Maf1_MB_CHTN_PM_2_depth$AF>=filter)
str(Maf1_MB_CHTN_PM_2_depth_AF3)

filter = 1
Maf1_MB_CHTN_PM_2_depth_AF1 = subset(Maf1_MB_CHTN_PM_2_depth,Maf1_MB_CHTN_PM_2_depth$AF>=filter)
str(Maf1_MB_CHTN_PM_2_depth_AF1)

filter = 0.5
Maf1_MB_CHTN_PM_2_depth_AF = subset(Maf1_MB_CHTN_PM_2_depth,Maf1_MB_CHTN_PM_2_depth$AF>=filter)
str(Maf1_MB_CHTN_PM_2_depth_AF)
cat('Maf1_MB_CHTN_PM_2_depth_AF in dbSNP: ',length(which(! Maf1_MB_CHTN_PM_2_depth_AF$dbSNP_RS=='novel')))
cat(' Maf1_MB_CHTN_PM_2_depth_AF NOT dbSNP: ',length(which(Maf1_MB_CHTN_PM_2_depth_AF$dbSNP_RS=='novel')))

Maf1_MB_CHTN_PM_6 <- subsetMaf(maf1, tsb = "ABSJ052018CHTN-6-PM"  , genes = NULL, fields = NULL, query = NULL,
                               mafObj = FALSE, includeSyn = TRUE, isTCGA = FALSE, restrictTo = "all")

min.depth.filter  = 100
Maf1_MB_CHTN_PM_6_depth = subset(Maf1_MB_CHTN_PM_6,Maf1_MB_CHTN_PM_6$t_depth>=min.depth.filter)
str(Maf1_MB_CHTN_PM_6_depth)

cat('Maf1_MB_CHTN_PM_6_depth in dbSNP: ',length(which(! Maf1_MB_CHTN_PM_6_depth$dbSNP_RS=='novel')))
cat('Maf1_MB_CHTN_PM_6_depth NOT dbSNP: ',length(which(Maf1_MB_CHTN_PM_6_depth$dbSNP_RS=='novel')))


filter = 3
Maf1_MB_CHTN_PM_6_depth_AF3 = subset(Maf1_MB_CHTN_PM_6_depth,Maf1_MB_CHTN_PM_6_depth$AF>=filter)
str(Maf1_MB_CHTN_PM_6_depth_AF3)

filter = 1
Maf1_MB_CHTN_PM_6_depth_AF1 = subset(Maf1_MB_CHTN_PM_6_depth,Maf1_MB_CHTN_PM_6_depth$AF>=filter)
str(Maf1_MB_CHTN_PM_6_depth_AF1)

filter = 0.5
Maf1_MB_CHTN_PM_6_depth_AF = subset(Maf1_MB_CHTN_PM_6_depth,Maf1_MB_CHTN_PM_6_depth$AF>=filter)
str(Maf1_MB_CHTN_PM_6_depth_AF)
cat('Maf1_MB_CHTN_PM_6_depth_AF in dbSNP: ',length(which(! Maf1_MB_CHTN_PM_6_depth_AF$dbSNP_RS=='novel')))
cat(' Maf1_MB_CHTN_PM_6_depth_AF NOT dbSNP: ',length(which(Maf1_MB_CHTN_PM_6_depth_AF$dbSNP_RS=='novel')))

Maf1_MB_CHTN_PM_4 <- subsetMaf(maf1, tsb = "ABSJ052018CHTN-4-PM"  , genes = NULL, fields = NULL, query = NULL,
                               mafObj = FALSE, includeSyn = TRUE, isTCGA = FALSE, restrictTo = "all")

min.depth.filter  = 100
Maf1_MB_CHTN_PM_4_depth = subset(Maf1_MB_CHTN_PM_4,Maf1_MB_CHTN_PM_4$t_depth>=min.depth.filter)
str(Maf1_MB_CHTN_PM_4_depth)

cat('Maf1_MB_CHTN_PM_4_depth in dbSNP: ',length(which(! Maf1_MB_CHTN_PM_4_depth$dbSNP_RS=='novel')))
cat('Maf1_MB_CHTN_PM_4_depth NOT dbSNP: ',length(which(Maf1_MB_CHTN_PM_4_depth$dbSNP_RS=='novel')))

filter = 3
Maf1_MB_CHTN_PM_4_depth_AF3 = subset(Maf1_MB_CHTN_PM_4_depth,Maf1_MB_CHTN_PM_4_depth$AF>=filter)
str(Maf1_MB_CHTN_PM_4_depth_AF3)

filter = 1
Maf1_MB_CHTN_PM_4_depth_AF1 = subset(Maf1_MB_CHTN_PM_4_depth,Maf1_MB_CHTN_PM_4_depth$AF>=filter)
str(Maf1_MB_CHTN_PM_4_depth_AF1)

filter = 0.5
Maf1_MB_CHTN_PM_4_depth_AF = subset(Maf1_MB_CHTN_PM_4_depth,Maf1_MB_CHTN_PM_4_depth$AF>=filter)
str(Maf1_MB_CHTN_PM_4_depth_AF)
cat('Maf1_MB_CHTN_PM_4_depth_AF in dbSNP: ',length(which(! Maf1_MB_CHTN_PM_4_depth_AF$dbSNP_RS=='novel')))
cat(' Maf1_MB_CHTN_PM_4_depth_AF NOT dbSNP: ',length(which(Maf1_MB_CHTN_PM_4_depth_AF$dbSNP_RS=='novel')))


Maf1_MB_CHTN_PM_7 <- subsetMaf(maf1, tsb = "ABSJ052018CHTN-7-PM"  , genes = NULL, fields = NULL, query = NULL,
                               mafObj = FALSE, includeSyn = TRUE, isTCGA = FALSE, restrictTo = "all")

min.depth.filter  = 100
Maf1_MB_CHTN_PM_7_depth = subset(Maf1_MB_CHTN_PM_7,Maf1_MB_CHTN_PM_7$t_depth>=min.depth.filter)
str(Maf1_MB_CHTN_PM_7_depth)
cat('Maf1_MB_CHTN_PM_7_depth in dbSNP: ',length(which(! Maf1_MB_CHTN_PM_7_depth$dbSNP_RS=='novel')))
cat('Maf1_MB_CHTN_PM_7_depth NOT dbSNP: ',length(which(Maf1_MB_CHTN_PM_7_depth$dbSNP_RS=='novel')))




Maf1_MB_CHTN_PM_1 <- subsetMaf(maf1, tsb = "ABSJ052018CHTN-1-PM"  , genes = NULL, fields = NULL, query = NULL,
                               mafObj = FALSE, includeSyn = TRUE, isTCGA = FALSE, restrictTo = "all")

min.depth.filter  = 100
Maf1_MB_CHTN_PM_1_depth = subset(Maf1_MB_CHTN_PM_1,Maf1_MB_CHTN_PM_1$t_depth>=min.depth.filter)
str(Maf1_MB_CHTN_PM_7_depth)
cat('Maf1_MB_CHTN_PM_1_depth in dbSNP: ',length(which(! Maf1_MB_CHTN_PM_1_depth$dbSNP_RS=='novel')))
cat('Maf1_MB_CHTN_PM_1_depth NOT dbSNP: ',length(which(Maf1_MB_CHTN_PM_1_depth$dbSNP_RS=='novel')))

filter = 3
Maf1_MB_CHTN_PM_1_depth_AF3 = subset(Maf1_MB_CHTN_PM_1_depth,Maf1_MB_CHTN_PM_1_depth$AF>=filter)
str(Maf1_MB_CHTN_PM_1_depth_AF3)

filter = 1
Maf1_MB_CHTN_PM_1_depth_AF1 = subset(Maf1_MB_CHTN_PM_1_depth,Maf1_MB_CHTN_PM_1_depth$AF>=filter)
str(Maf1_MB_CHTN_PM_1_depth_AF1)

filter = 0.5
Maf1_MB_CHTN_PM_1_depth_AF = subset(Maf1_MB_CHTN_PM_1_depth,Maf1_MB_CHTN_PM_1_depth$AF>=filter)
str(Maf1_MB_CHTN_PM_1_depth_AF)
cat('Maf1_MB_CHTN_PM_1_depth_AF in dbSNP: ',length(which(! Maf1_MB_CHTN_PM_1_depth_AF$dbSNP_RS=='novel')))
cat(' Maf1_MB_CHTN_PM_1_depth_AF NOT dbSNP: ',length(which(Maf1_MB_CHTN_PM_1_depth_AF$dbSNP_RS=='novel')))

str(Maf1_df)
##658341 obs. of  132 variables
View(Maf1_df)
mean_cl_normal(Maf1_df$t_depth)
median_hilow(Maf1_df$t_depth)
boxplot(Maf1_df$dbSNP_RS ,Mafdf$t_depth)
hist(Maf1_df$t_depth)
hist(Maf1_df$t_depth,include.lowest = FALSE, freq = TRUE)
qplot(Maf1_df$t_depth)
min.depth.filter  = 30

ggplot(Maf1_df, aes(x = Maf1_df$t_depth)) + geom_histogram(aes(y = ..density..), binwidth = 20)

min.depth.filter  = 300
Maf1_df_depth = subset(Maf1_df,Maf1_df$t_depth>=min.depth.filter)
ggplot(Maf1_df_depth, aes(x = Maf1_df_depth$t_depth)) + geom_histogram(aes(y = ..count..), binwidth = 75)   

p <- ggplot(Maf1_df_depth, aes(x = Maf1_df_depth$AF,  y = Maf1_df_depth$t_depth))
p + geom_point()
## filtering for depth

min.depth.filter  = 150
Maf1_df_depth = subset(Maf1_df,Maf1_df$t_depth>=min.depth.filter)
str(Maf1_df_depth)

## Filetring for depth lead to  22433 obs. of  132 variables

min.depth.filter  = 300
Maf1_df_depth = subset(Maf1_df,Maf1_df$t_depth>=min.depth.filter)
str(Maf1_df_depth)

## Filetring for depth lead to  10554 obs. of  132 variables
min.depth.filter  = 30
Maf1_df_depth = subset(Maf1_df,Maf1_df$t_depth>=min.depth.filter)
str(Maf1_df_depth)


## calculating the percentage of variants that are prsent in dbsnp; higher ratio indicates contamination with germline

cat('Maf1_df_depth in dbSNP: ',length(which(! Maf1_df_depth$dbSNP_RS=='novel')))

##9419 (89.25%) of the variants are prsent in dbsnp

VAF_filter = 0.5
Maf1_df_depth_VAF = subset(Maf1_df_depth,Maf1_df_depth$gnomAD_AF>=filter)
str(Maf1_df_depth_VAF)

##	740 obs. of  132 variables

VAF_filter = 3
Maf1_df_depth_VAF = subset(Maf1_df_depth,Maf1_df_depth$gnomAD_AF>=filter)
str(Maf1_df_depth_VAF)

VAF_filter = 3
Maf1_df_depth_VAF = subset(Maf1_df_depth,Maf1_df_depth$AF>=filter)
str(Maf1_df_depth_VAF)

VAF_filter = 3
Maf1_df_depth_VAF_AF = subset(Maf1_df_depth_VAF,Maf1_df_depth_VAF$gnomAD_AF>=filter)
str(Maf1_df_depth_VAF_AF)
View(Maf1_df_depth_VAF_AF)




##


Maf1_MB2 <- subsetMaf(maf1, tsb = "ABSJ112017MB2", genes = NULL, fields = NULL, query = NULL,
                      mafObj = FALSE, includeSyn = TRUE, isTCGA = FALSE, restrictTo = "all")
str(Maf1_MB2)
head(Maf1_MB2)
data.frame(Maf1_MB2$gnomAD_AF,Maf1_MB2$Start_Position)
hist(Maf1_MB2$gnomAD_AF,breaks = 3,main = Maf1_MB2)



cat('Maf1_MB2 in dbSNP: ',length(which(! Maf1_MB2$dbSNP_RS=='novel')))
cat('Maf1_MB2 NOT dbSNP: ',length(which(Maf1_MB2$dbSNP_RS=='novel')))

min.depth.filter  = 
  
  Maf1_MB2_filter = subset(Maf1_MB2,Maf1_MB2$t_depth>=min.depth.filter)
str(Maf1_MB2_filter)

cat('Maf1_MB2_filter NOT dbSNP: ',length(which(Maf1_MB2_filter$dbSNP_RS=='novel')))

filter = 0.5
Maf1_MB2_filter_AF = subset(Maf1_MB2_filter,Maf1_MB2_filter$gnomAD_AF>=filter)
str(Maf1_MB2_filter_AF)
data.frame(Maf1_MB2_filter_AF$gnomAD_AF,Maf1_MB2_filter_AF$Start_Position)
p <- ggplot(Depth_dbsnp_df, aes(x = Depth_dbsnp_df$Tumor_Sample_Barcode,  y = Depth_dbsnp_df$Variants))
p + geom_point()

ggplot(Depth_dbsnp_df,aes(x = Depth_dbsnp_df$Tumor_Sample_Barcode, y = Depth_dbsnp_df$Novel)) + geom_col(fill = "black") 