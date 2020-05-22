# HZ19 qPCRs

library(Rmisc)
library(httr)
library(RCurl)
library(plyr)
library(dplyr)
library(reshape2)
library(tidyverse)

qPCR1 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/raw_qPCR/HZ19_qPCR/HZ19_qPCR1.csv"))
qPCR2 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/raw_qPCR/HZ19_qPCR/HZ19_qPCR2.csv"))
qPCR3 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/raw_qPCR/HZ19_qPCR/HZ19_qPCR3.csv"))
qPCR4 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/raw_qPCR/HZ19_qPCR/HZ19_qPCR4.csv"))
qPCR5 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/raw_qPCR/HZ19_qPCR/HZ19_qPCR5.csv"))
qPCR6 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/raw_qPCR/HZ19_qPCR/HZ19_qPCR6.csv"))

qPCR <- rbind(qPCR1, qPCR2)
qPCR <- rbind(qPCR, qPCR3)
qPCR <- rbind(qPCR, qPCR4)
qPCR <- rbind(qPCR, qPCR5)
qPCR <- rbind(qPCR, qPCR6)
qPCR$Mouse_ID <- as.character(qPCR$Mouse_ID)

names(qPCR)[names(qPCR) == "Detector"] <- "Target"
names(qPCR)[names(qPCR) == "Sample.Name"] <- "Mouse_ID"
qPCR$Ct <- as.numeric(as.character(qPCR$Ct))

qPCR.long <- qPCR %>% dplyr::group_by(Mouse_ID, Target) %>% dplyr::summarise(Ct = mean(Ct, na.rm = T))
qPCR.long <- data.frame(qPCR.long)
qPCR.wide <- reshape(qPCR.long[, c("Target", "Mouse_ID","Ct")],
                   timevar = "Target", idvar = "Mouse_ID", direction = "wide")

qPCR.wide$delta <- qPCR.wide$Ct.Mus - qPCR.wide$Ct.Eimeria
qPCR.wide$Mouse_ID <- sub("^", "AA_0", qPCR.wide$Mouse_ID)
qPCR.wide$EXP_type <- "wild"

write.csv(qPCR.wide, "/Users/Luke Bednar/Mouse_Eimeria_Databasing/data/Eimeria_detection/HZ19_qPCR.csv")
