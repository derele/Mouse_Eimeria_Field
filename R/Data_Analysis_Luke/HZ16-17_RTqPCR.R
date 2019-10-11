library(httr)
library(RCurl)
library(tidyverse)
library(Rmisc)

#read in 2016-2017 qPCRs
PCR <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/FINALqpcrData_2016_2017_threshold5.csv"
PCR <- read.csv(text = getURL(PCR))
# select positive samples
positive2016 <- filter(PCR, qPCRstatus == "positive", year == 2016, qPCRsummary == "infected cecum")
positive2017 <- filter(PCR, qPCRstatus == "positive", year == 2017, qPCRsummary == "infected cecum")
# select negative samples
negative2016 <- filter(PCR, qPCRstatus == "negative", year == 2016, qPCRsummary == "non infected")
negative2017 <- filter(PCR, qPCRstatus == "negative", year == 2017, qPCRsummary == "non infected")
# randomly select negative samples (21 per 2016, 29 per 2017)
negative2016 <- sample_n(negative2016, 21)
negative2017 <- sample_n(negative2017, 29)
#merge and save
negative <- bind_rows(negative2016, negative2017)
positive <- bind_rows(positive2016, positive2017)
HZ16and17 <- bind_rows(negative, positive)
write.csv(HZ16and17, "//home/luke/Repositories/Mouse_Eimeria_Databasing/Mouse_Eimeria_Databasing/data/HZ16and17RT-qPCRsample_list.csv")
