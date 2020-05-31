library(Rmisc)
library(httr)
library(RCurl)
library(plyr)
library(dplyr)
library(reshape2)
library(tidyverse)

qPCR <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/Svenja/table_ct_and_more.csv"))
names(qPCR)[names(qPCR) == "Name"] <- "Mouse_ID"
qPCR <- qPCR %>% separate(Mouse_ID, c("CEWE", "AA", "Mouse_ID"))
qPCR$Mouse_ID <- sub("^", "AA_0", qPCR$Mouse_ID )
qPCR$CEWE <- NULL
qPCR$AA <- NULL
qPCR <- select(qPCR, Mouse_ID, Eimeria.presence.in.Caecum, deltaCtMmE_tissue)
names(qPCR)[names(qPCR) == "Eimeria.presence.in.Caecum"] <- "MC.Eimeria"
names(qPCR)[names(qPCR) == "deltaCtMmE_tissue"] <- "delta"

write.csv(qPCR, "/Users/Luke Bednar/Mouse_Eimeria_Databasing/data/Eimeria_detection/HZ18_qPCR.csv")
