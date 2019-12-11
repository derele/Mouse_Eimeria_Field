library(httr)
library(RCurl)
library(Rmisc)
library(dplyr)

#read in 2016-2017 qPCRs (based on threshold, check melting curves too)
Mice <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/MiceTable_fullEimeriaInfos_2014to2017.csv"
Mice <- read.csv(text = getURL(Mice))
# reduce this monstrosity
Mice <- dplyr::select(Mice, Mouse_ID, Sex, Longitude, Latitude, Year, HI_NLoci, HI, Body_weight, Body_length, Species,
                      Spleen_mass, Status, mean_neubauer, OPG, delta_ct_cewe_MminusE, eimeriaSpecies, observer)
E64 <- dplyr::filter(Mice, eimeriaSpecies == "E_ferrisi")
E88 <- dplyr::filter(Mice, eimeriaSpecies == "E_falciformis")
Double <- dplyr::filter(Mice, eimeriaSpecies == "Double")
Double_tbd <- dplyr::filter(Mice, eimeriaSpecies == "Double_tbd")
Double_det <- dplyr::filter(Mice, eimeriaSpecies == "Double_ferrisi_vermiformis")
Other <- dplyr::filter(Mice, eimeriaSpecies == "Other")
Negative <- dplyr::filter(Mice, eimeriaSpecies == "Negative")

Positive <- rbind(E64, E88)
Positive <- rbind(Positive, Double)
Positive <- rbind(Positive, Double_tbd)
Positive <- rbind(Positive, Double_det)
Positive <- rbind(Positive, Other)

write.csv(Positive, "~/Repositories/Mouse_Eimeria_Databasing/Mouse_Eimeria_Databasing/data/Eimeria_detection/HZ14-17_positives.csv")
write.csv(Negative, "~/Repositories/Mouse_Eimeria_Databasing/Mouse_Eimeria_Databasing/data/Eimeria_detection/HZ14-17_negatives.csv")







# 
# 
# cDNA <- read.csv("./Downloads/cDNA wild mice 16_17 new.CSV")
# cDNA <- select(cDNA, Sample.ID, Nucleic.Acid.Conc.)
# colnames(cDNA)[1] <- "Mouse_ID"
# 
# positives <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/MC_verified_positives.csv"
# positives <- read.csv(text = getURL(positives))
# positives <- select(positives, Mouse_ID)
# 
# cDNA <- merge(cDNA, positives, by = "Mouse_ID")
# 
# write.csv(cDNA, file = "./Desktop/HZ16-17_cDNA_concentrations.csv")
# 
# HZ16_17 <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/HZ16-17_InfInt_MC_Lorenzo%26Mert.csv"
# HZ16_17 <- read.csv(text = getURL(HZ16_17))
# HZ16_17 <- select(HZ16_17, Mouse_ID, observer)
# HZ16_17 <- filter(HZ16_17, observer == "Lorenzo")
# HZ16_17 <- merge(HZ16_17, cDNA)
# HZ16_17 - positives
