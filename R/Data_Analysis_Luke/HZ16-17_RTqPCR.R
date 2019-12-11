library(httr)
library(RCurl)
library(Rmisc)
library(tidyverse)
library(purrr)
library(ggplot2)

# load in sample list


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

