library(Rmisc)
library(httr)
library(RCurl)
library(plyr)
library(dplyr)
library(reshape2)
library(tidyverse)

qPCR <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/HZ19_qPCR.csv"))
qPCR$X <- NULL
# doesn"t exist yet
#RT <- read.csv(text = getURL(""))

ELISA_CEWE <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/ELISAs/HZ19_CEWE_ELISAs_complete.csv"))
ELISA_CEWE$X <- NULL

FACS <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Field_data/HZ19_FACS_complete.csv"))
FACS$X <- NULL

immuno <- merge(qPCR, ELISA_CEWE, all = T)
immuno <- merge(immuno, FACS, all = T)

write.csv(immuno, "/Users/Luke Bednar/Mouse_Eimeria_Databasing/data/HZ19_immuno.csv")
