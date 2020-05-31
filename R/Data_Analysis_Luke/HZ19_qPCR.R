# HZ19 qPCRs

library(Rmisc)
library(httr)
library(RCurl)
library(plyr)
library(dplyr)
library(reshape2)
library(tidyverse)
# load in raw tables
qPCR1 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/raw_qPCR/HZ19_qPCR/HZ19_qPCR1_clean.csv"))
qPCR2 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/raw_qPCR/HZ19_qPCR/HZ19_qPCR2_clean.csv"))
qPCR3 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/raw_qPCR/HZ19_qPCR/HZ19_qPCR3_clean.csv"))
qPCR4 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/raw_qPCR/HZ19_qPCR/HZ19_qPCR4_clean.csv"))
qPCR5 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/raw_qPCR/HZ19_qPCR/HZ19_qPCR5_clean.csv"))
qPCR6 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/raw_qPCR/HZ19_qPCR/HZ19_qPCR6_clean.csv"))
qPCR7 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/raw_qPCR/HZ19_qPCR/HZ19_qPCR7_clean.csv"))
qPCR8 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/raw_qPCR/HZ19_qPCR/HZ19_qPCR8_clean.csv"))
qPCR9 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/raw_qPCR/HZ19_qPCR/HZ19_qPCR9_clean.csv"))
qPCR10 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/raw_qPCR/HZ19_qPCR/HZ19_qPCR10_clean.csv"))
qPCR11 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/raw_qPCR/HZ19_qPCR/HZ19_qPCR11_clean.csv"))
# easier to add this and then remove Well column from all once rbound
qPCR1$Well <- "a"

qPCR <- rbind(qPCR1, qPCR2)
qPCR <- rbind(qPCR, qPCR3)
qPCR <- rbind(qPCR, qPCR4)
qPCR <- rbind(qPCR, qPCR5)
qPCR <- rbind(qPCR, qPCR6)
qPCR <- rbind(qPCR, qPCR7)
qPCR <- rbind(qPCR, qPCR8)
qPCR <- rbind(qPCR, qPCR9)
qPCR <- rbind(qPCR, qPCR10)
qPCR <- rbind(qPCR, qPCR11)

# set names and convert columns to appropriate type
names(qPCR)[names(qPCR) == "Detector"] <- "Target"
names(qPCR)[names(qPCR) == "Sample.Name"] <- "Mouse_ID"
qPCR$Ct <- as.numeric(as.character(qPCR$Ct))
qPCR$Mouse_ID <- as.character(qPCR$Mouse_ID)
# remove well numebrs and NTCs
qPCR$Well <- NULL
qPCR <- subset(qPCR, !qPCR$Mouse_ID == "NTC")


qPCR.long <- qPCR %>% dplyr::group_by(Mouse_ID, Target, MC) %>% dplyr::summarise(Ct = mean(Ct, na.rm = T))
qPCR.long <- data.frame(qPCR.long)
qPCR.wide <- reshape(qPCR.long[, c("Target", "Mouse_ID","Ct", "MC")],
                   timevar = "Target", idvar = "Mouse_ID", direction = "wide")

qPCR.wide$delta <- qPCR.wide$Ct.Mus - qPCR.wide$Ct.Eimeria
qPCR.wide$Mouse_ID <- sub("^", "AA_0", qPCR.wide$Mouse_ID)
qPCR.wide$EXP_type <- "wild"
# remove Mus MC because they are all true
qPCR.wide$MC.Mus <- NULL

write.csv(qPCR.wide, "/Users/Luke Bednar/Mouse_Eimeria_Databasing/data/Eimeria_detection/HZ19_qPCR.csv")
