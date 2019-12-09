library(httr)
library(RCurl)
library(Rmisc)
library(tidyverse)
library(purrr)
library(ggplot2)

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

write.csv(Positive, "~/Mouse_Eimeria_Databasing/data/Eimeria_detection/HZ14-17_positives.csv")

### Check against MC analysis (Lorenzo)
LorenzoMC <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/HZ16-17_InfInt_MC_Lorenzo%26Mert.csv"
LorenzoMC <- read.csv(text = getURL(LorenzoMC))

LorenzoMC <- merge(LorenzoMC, Positive)
TruePositives <- subset(LorenzoMC, Caecum == "pos")
write.csv(TruePositives, file = "~/Mouse_Eimeria_Databasing/Mouse_Eimeria_Databasing/data/Eimeria_detection/MC_verified_positives.csv")

############################## Add oocyst data
oocysts <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/Eimeria_oocysts_2015%262017_Lorenzo.csv"
oocysts <- read.csv(text = getURL(oocysts))
oocysts <- select(oocysts, Mouse_ID, mean_neubauer, OPG)
HZ16and17 <- merge(oocysts, HZ16and17, by = "Mouse_ID", all.y = TRUE)
# check for oocyst positive and qPCR negative flukes (keeps reducing the overall number... maybe stick to oocyst and HZ merge)
DoublePositive <- dplyr::filter(HZ16and17, qPCRstatus == "positive", OPG > 0)
DoublePositive$positive <- "double"
OocystPositive <- dplyr::filter(HZ16and17, qPCRstatus == "negative", OPG > 0)
OocystPositive$positive <- "oocyst"
qPCRPositive <- dplyr::filter(HZ16and17, qPCRstatus == "positive", OPG == 0)
qPCRPositive$positive <- "qPCR"
PositiveInvestigate <- rbind(DoublePositive, OocystPositive, qPCRPositive)
PositiveInvestigate <- distinct(PositiveInvestigate)
HZ16and17 <- merge(HZ16and17, PositiveInvestigate, by = "Mouse_ID")

################### Add melting cuve analysis data #######################

