library(httr)
library(RCurl)
library(tidyverse)
library(Rmisc)
library(purrr)

#read in 2016-2017 qPCRs (based on threshold, check emlting curves too)
Mice <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/MiceTable_fullEimeriaInfos_2014to2017.csv"
Mice <- read.csv(text = getURL(Mice))
# reduce this monstrosity
Mice <- dplyr::select(Mice, Mouse_ID, Sex, Longitude, Latitude, Year, HI_NLoci, HI, Body_weight, Body_length, Species,
                      Spleen_mass, Status, mean_neubauer, OPG, delta_ct_cewe_MminusE, eimeriaSpecies)
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



data.frame(E64, E88, Double, Double_tbd, Double_det, Other) %>% reduce(inner_join, by = "Mouse_ID")

# select positive samples 2016
positive2016a <- filter(PCR, qPCRstatus == "positive", year == 2016, qPCRsummary == "infected cecum")
positive2016b <- filter(PCR, qPCRstatus == "positive", year == 2016, qPCRsummary == "cecum stronger")
positive2016c <- filter(PCR, qPCRstatus == "positive", year == 2016, qPCRsummary == "ileum stronger")
# select positive samples 2017
positive2017a <- filter(PCR, qPCRstatus == "positive", year == 2017, qPCRsummary == "infected cecum")
positive2017b <- filter(PCR, qPCRstatus == "positive", year == 2017, qPCRsummary == "cecum stronger")
positive2017c <- filter(PCR, qPCRstatus == "positive", year == 2017, qPCRsummary == "ileum stronger")
# select negative samples
negative2016 <- filter(PCR, qPCRstatus == "negative", year == 2016, qPCRsummary == "non infected")
negative2017 <- filter(PCR, qPCRstatus == "negative", year == 2017, qPCRsummary == "non infected")
# randomly select negative samples (28 per 2016, 29 per 2017)
negative2016 <- sample_n(negative2016, 28)
negative2017 <- sample_n(negative2017, 33)
negative <- bind_rows(negative2016, negative2017)
#merge
positiveA <- bind_rows(positive2016a, positive2017a)
positiveB <- bind_rows(positive2016b, positive2017b)
positiveC <- bind_rows(positive2016c, positive2017c)
positiveAB <- bind_rows(positiveA, positiveB)
positive <- bind_rows(positiveAB, positiveC)
dplyr::arrange(positive, year)
HZ16and17 <- rbind(positive, negative)
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

