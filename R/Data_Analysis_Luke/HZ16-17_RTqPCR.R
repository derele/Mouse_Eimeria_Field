library(httr)
library(RCurl)
library(tidyverse)
library(Rmisc)

#read in 2016-2017 qPCRs (based on threshold, check emlting curves too)
PCR <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/FINALqpcrData_2016_2017_threshold5.csv"
PCR <- read.csv(text = getURL(PCR))
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

