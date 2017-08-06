#########################
library(ggplot2)
library(gridExtra)
library(reshape2)
library(scales)
library(RCurl)
library(ggrepel)

#########################
## Import:
rawdata2016 <- read.csv("../data_raw/Results_flotation_and_PCR_2016.csv")

## Remove AA_0001 to AA_0046 that was 1 farm done before the trip :
data2016 <- rawdata2016[-(1:46),]

## Observation : who performed the counting
table(data2016$Observation)

#########################
# Selection : data POSIIVE for flotation by standard methods +
# negative for flotation & positive for AP5 & positive for all other 3 markers :
  # data NEGATIVE if negative for flotation by standard methods & negative for AP5 +
# negative for flotation & positive for AP5 & negative for all other 3 markers :

# Function to define the % positive :
prevalence <- function(data){
  positive <- sum(na.omit(data$Flotation == "POSITIVE")) +
    sum(na.omit(data$Flotation == "NEGATIVE" & data$PCR_AP5 == "POSITIVE" & 
                  data$at.least.1.other.marker_ORF_COI_or_18S == "POSITIVE"))
  
  total <- sum(!is.na(data$Flotation))
  
  positive/total*100
}

## What are the discrepencies?
# 1. standard methods of flotation :
prevalence(data2016[which(data2016$Observation == "P"),])

# 2. Non standard : same but do not consider also non standard method (to compare)
prevalence(data2016)

# 3. Considering only samples still there :
prevalence(data2016[which(!is.na(data2016$Box_No)),])

# 4. Detection by AP5 : (to compare)
table(data2016$PCR_AP5)[names(table(data2016$PCR_AP5)) == "POSITIVE"]/
  sum(table(data2016$PCR_AP5))*100

#########################
## Improve our data :

# BEWARE sample AA_0094 counted twice

## A. needs to recount :
#write.csv(x = data.frame(sample = rawdata2016[which(rawdata2016$A.need.recount),]$Sample_ID,
#                         box = rawdata2016[which(rawdata2016$A.need.recount),]$Box_No), file = "alice.csv")

## E. needs to find boxes :
#write.csv(x = data.frame(missing.samples.2016 = rawdata2016[which(rawdata2016$find_sample),]$Sample_ID), file = "missing.csv")

## P. needs to do some AP5 :
#write.csv(x = data.frame(AP5_to_do = rawdata2016[which(rawdata2016$AP5_to_do),]$Sample_ID), file = "phuong.csv")

## V. needs to do some PCR other markers :
#write.csv(x = data.frame(othermarkers_to_do = rawdata2016[which(rawdata2016$other_markers_to_test),]$Sample_ID),
#                         file = "victor.csv")

## Plot BoxOrNoBox + observer :
data2016$BoxOrNoBox <- "FALSE"
data2016$BoxOrNoBox[which(!is.na(data2016$Box_No))] <- "TRUE"

data2016$Observer_and_box <- paste(data2016$Observation, data2016$BoxOrNoBox)

ggplot(data = data2016 ,
       aes(x = as.factor(Sample_ID), y = OPG_.oocysts.ml.g_of_faeces.,
           col = Observer_and_box)) +
  scale_color_manual(values = c("red", "green", "orange", "blue")) +
  geom_point(size = 3)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#########################
## Export :
# 1. check that we have flotation for all positive samples :
dataexport <- data2016[which(data2016$Observation == "P"),]

prevalence(dataexport)

table(dataexport$OPG_.oocysts.ml.g_of_faeces. != 0)[which(names(table(dataexport$OPG_.oocysts.ml.g_of_faeces. != 0)) == "TRUE")]/
  sum(table(dataexport$OPG_.oocysts.ml.g_of_faeces.)) * 100

#2. If not the same, rely on the counts will under-estimate the positive (no other choice though)
dataexport <- data.frame(Mouse_ID = dataexport$Sample_ID,
                         OPG = dataexport$OPG_.oocysts.ml.g_of_faeces.)

write.csv(x = dataexport, file = "../data_clean/Results_flotation_and_PCR_2016_CLEAN.csv", row.names = FALSE)
