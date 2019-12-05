library(httr)
library(RCurl)
library(Rmisc)
library(dplyr)
cDNA <- read.csv("./Downloads/cDNA wild mice 16_17 new.CSV")
cDNA <- select(cDNA, Sample.ID, Nucleic.Acid.Conc.)
colnames(cDNA)[1] <- "Mouse_ID"

positives <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/MC_verified_positives.csv"
positives <- read.csv(text = getURL(positives))
positives <- select(positives, Mouse_ID)

cDNA <- merge(cDNA, positives, by = "Mouse_ID")

write.csv(cDNA, file = "./Desktop/HZ16-17_cDNA_concentrations.csv")

HZ16_17 <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/HZ16-17_InfInt_MC_Lorenzo%26Mert.csv"
HZ16_17 <- read.csv(text = getURL(HZ16_17))
HZ16_17 <- select(HZ16_17, Mouse_ID, observer)
HZ16_17 <- filter(HZ16_17, observer == "Lorenzo")
HZ16_17 <- merge(HZ16_17, cDNA)
HZ16_17 - positives
