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

## Oocyst per gram calculation :
data2016$Raw_count_by_A_Neubauer_1bigsquare <- as.numeric(as.character(data2016$Raw_count_by_A_Neubauer_1bigsquare))
data2016$Ooperg.A <- (data2016$Raw_count_by_A_Neubauer_1bigsquare * 10000 * data2016$Volume_mL)/data2016$Weight_feces_.g.

data2016$Raw_count_by_P_Neubauer_1bigsquare <- as.numeric(as.character(data2016$Raw_count_by_P_Neubauer_1bigsquare))
data2016$Ooperg.P <- (data2016$Raw_count_by_P_Neubauer_1bigsquare * 10000 * data2016$Volume_mL)/data2016$Weight_feces_.g.

## Compare the raw counts between 3 measurers :
ggplot(data = data2016, aes(x = Sample_ID)) +
  geom_point(aes(y = data2016$Raw_count_by_A_Neubauer_1bigsquare), col = "pink") + 
  geom_point(aes(y = data2016$Raw_count_by_P_Neubauer_1bigsquare), col = "blue") + 
  geom_point(aes(y = data2016$Raw_count_by_E), col = "orange") 


(data2016$Raw_count_by_E * data2016$Volume_mL) / data2016$Raw_count_by_A_Neubauer_1bigsquare 

data2016$Raw_count_by_P_Neubauer_1bigsquare / data2016$Raw_count_by_A_Neubauer_1bigsquare 

# samples to chack by other markers :
data2016$PCR_Results_Ap5 != data2016$Raw_count_by_A_Neubauer_1bigsquare

# how many mice?
length(unique(data2016$Sample_ID))

# how many mice that we still have the sample?
datawegot <- data2016[which(!is.na(data2016$Box_No)),]

# how many samples were negative for Ap5 and flotation and where discarded?
datanull <- data2016[which(is.na(data2016$Box_No)),]
datanull <- datanull[which(datanull$Flotation == "NEGATIVE" & datanull$PCR_AP5 == "NEGATIVE"), ]
length(unique(datanull$Sample_ID))

# these are the samples to be used in the study :
data.final.2016 <- rbind(datawegot, datanull)
length(unique(data.final.2016$Sample_ID))

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
prevalence(datawegot)

# 4. Detection by AP5 : (to compare)
table(data2016$PCR_AP5)[names(table(data2016$PCR_AP5)) == "POSITIVE"]/
  sum(table(data2016$PCR_AP5))*100

# 5. the sample we will use :
prevalence(data.final.2016)

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

## Prevalence if we consider original data :
names(data2016)[1]<- "Mouse_ID"
ex <- merge(data2016, getGenDF(), by = "Mouse_ID")
data.agg <- aggregate(ex$OPG_.oocysts.ml.g_of_faeces., by = list(ex$Code), 
                      FUN = sum, na.rm=TRUE)

length(data.agg[which(data.agg$x != 0), ]$Group.1) / length(data.agg$Group.1) *100

#########################
## Export :
dataexport <- data.frame(Mouse_ID = data.final.2016$Sample_ID,
                         OPG = data.final.2016$OPG_.oocysts.ml.g_of_faeces.,
                         Year = 2016)

write.csv(x = dataexport, file = "../data_clean/Results_flotation_and_PCR_2016_CLEAN.csv", row.names = FALSE)