###### 2016 ######
qpcrData <- read.csv("../qPCR_2016.csv")

# Did Enas calculate the other way around?
qpcrData$delta_ct_cewe[qpcrData$observer_qpcr == "Enas"] <-
  - qpcrData$delta_ct_cewe[qpcrData$observer_qpcr == "Enas"]

qpcrData$delta_ct_ilwe[qpcrData$observer_qpcr == "Enas"] <-
  - qpcrData$delta_ct_ilwe[qpcrData$observer_qpcr == "Enas"]

# Here deltaCT = ct eimeria - ct mouse. If high infection, low deltaCT
# -deltaCT = ct mouse - ct eimeria
qpcrData$qPCRsummary[qpcrData$delta_ct_cewe > 6 & qpcrData$delta_ct_ilwe > 6] <- "non infected"
qpcrData$qPCRsummary[qpcrData$delta_ct_cewe < 6 & qpcrData$delta_ct_ilwe > 6] <- "infected cecum"
qpcrData$qPCRsummary[qpcrData$delta_ct_cewe > 6 & qpcrData$delta_ct_ilwe < 6] <- "infected ileum"

qpcrData$qPCRsummary[
  qpcrData$delta_ct_cewe < 6 & 
    qpcrData$delta_ct_ilwe < 6 & 
    qpcrData$delta_ct_cewe < qpcrData$delta_ct_ilwe] <- "cecum stronger"
qpcrData$qPCRsummary[
  qpcrData$delta_ct_cewe < 6 & 
    qpcrData$delta_ct_ilwe < 6 & 
    qpcrData$delta_ct_cewe > qpcrData$delta_ct_ilwe] <- "ileum stronger"

# Infected or not?
qpcrData$qPCRstatus <- "positive"
qpcrData$qPCRstatus[is.na(qpcrData$qPCRsummary)] <- NA
qpcrData$qPCRstatus[qpcrData$qPCRsummary %in% "non infected"] <- "negative"

# and keep the infected segment value OR the higher value 
qpcrData$delta_ct[
  qpcrData$qPCRsummary %in% c("infected cecum", "cecum stronger")] <- 
  qpcrData$delta_ct_cewe[
    qpcrData$qPCRsummary %in% c("infected cecum", "cecum stronger")] 

qpcrData$delta_ct[
  qpcrData$qPCRsummary %in% c("infected ileum", "ileum stronger")] <- 
  qpcrData$delta_ct_ilwe[
    qpcrData$qPCRsummary %in% c("infected ileum", "ileum stronger")] 

# Turn around
qpcrData$delta_ct_MminusE <- - qpcrData$delta_ct

# Set floor values
qpcrData$delta_ct_MminusE[is.na(qpcrData$delta_ct_MminusE)] <- -6

## Add Tabea's values!
rawData <- read.csv("TabeaRaw/Eimeria_qPCR_Tissue_220818.csv")

# Correct name
rawData[rawData$Name %in% "ILWE_AA_140", "Name"] <- "ILWE_AA_0140"

rawMeltData <- read.csv("TabeaRaw/Eimeria_qPCR_Tissue_220818_Melting.csv")

# Merge by Name and Pos
mergedDF <- merge(rawData, rawMeltData, by = c("Name", "Pos"))

# Remove useless lines
mergedDF <- mergedDF[!is.na(mergedDF$Ct.Mean.SYBR) & !mergedDF$Name %in% "NTC_Mouse", ]

# Add tissue and Mouse_ID
x <- strsplit(as.character(mergedDF$Name), "_", 1)
mergedDF$tissue <- sapply( x, "[", 1)
mergedDF$Mouse_ID <- paste0("AA_", sapply( x, "[", 3))
rm(x)

# calculate deltaCtMminusE
calculateDeltaCt <- function(df){
  sumDataMouse <- df[df$Target.SYBR %in% "mouse",]
  sumDataEimeria <- df[df$Target.SYBR %in% "eimeria",]
  mergedData <- merge(sumDataEimeria, sumDataMouse,
                      by = c("Mouse_ID", "tissue", "Name"))
  mergedData <- unique(mergedData)
  mergedData$deltaCt_MminusE <- as.numeric(as.character(mergedData$Ct.Mean.SYBR.y)) - 
    as.numeric(as.character(mergedData$Ct.Mean.SYBR.x)) # DeltaCt MOUSE minus EIMERIA
  return(mergedData)
}

mergedDF <- calculateDeltaCt(mergedDF)

mergedDF <- unique(mergedDF[c("Mouse_ID", "tissue", "deltaCt_MminusE")])

mergedDF$year <- 2016
mergedDF$observer_qPCR <- "Tabea"

finalDataClean <- mergedDF["Mouse_ID"]

# Add CEWE
finalDataClean <- merge(finalDataClean,
                        mergedDF[mergedDF$tissue %in% "CEWE", c("Mouse_ID", "deltaCt_MminusE")],
                        all.x = T)
names(finalDataClean)[names(finalDataClean) %in% "deltaCt_MminusE"] <- "delta_ct_cewe_MminusE"

# Add ILWE
finalDataClean <- merge(finalDataClean,
                        finalData[finalData$tissue %in% "ILWE", c("Mouse_ID", "deltaCt_MminusE")],
                        all.x = T)
names(finalDataClean)[names(finalDataClean) %in% "deltaCt_MminusE"] <- "delta_ct_ilwe_MminusE"

finalDataClean <- unique(finalDataClean)

# Add observer
finalDataClean$observer_qpcr <- "Lorenzo"

# Add year
finalDataClean$Year <- 2017




## TO TURN AROUND!!
# Here deltaCT = ct eimeria - ct mouse. If high infection, low deltaCT
# -deltaCT = ct mouse - ct eimeria
qpcrData$qPCRsummary[qpcrData$delta_ct_cewe > 6 & qpcrData$delta_ct_ilwe > 6] <- "non infected"
qpcrData$qPCRsummary[qpcrData$delta_ct_cewe < 6 & qpcrData$delta_ct_ilwe > 6] <- "infected cecum"
qpcrData$qPCRsummary[qpcrData$delta_ct_cewe > 6 & qpcrData$delta_ct_ilwe < 6] <- "infected ileum"

qpcrData$qPCRsummary[
  qpcrData$delta_ct_cewe < 6 & 
    qpcrData$delta_ct_ilwe < 6 & 
    qpcrData$delta_ct_cewe < qpcrData$delta_ct_ilwe] <- "cecum stronger"
qpcrData$qPCRsummary[
  qpcrData$delta_ct_cewe < 6 & 
    qpcrData$delta_ct_ilwe < 6 & 
    qpcrData$delta_ct_cewe > qpcrData$delta_ct_ilwe] <- "ileum stronger"

# Infected or not?
qpcrData$qPCRstatus <- "positive"
qpcrData$qPCRstatus[is.na(qpcrData$qPCRsummary)] <- NA
qpcrData$qPCRstatus[qpcrData$qPCRsummary %in% "non infected"] <- "negative"

# and keep the infected segment value OR the higher value 
qpcrData$delta_ct[
  qpcrData$qPCRsummary %in% c("infected cecum", "cecum stronger")] <- 
  qpcrData$delta_ct_cewe[
    qpcrData$qPCRsummary %in% c("infected cecum", "cecum stronger")] 

qpcrData$delta_ct[
  qpcrData$qPCRsummary %in% c("infected ileum", "ileum stronger")] <- 
  qpcrData$delta_ct_ilwe[
    qpcrData$qPCRsummary %in% c("infected ileum", "ileum stronger")] 






###### 2017 ######

# import data
rawData <- read.csv("./LorenzoRAW/CSVFiles/TotalLorenzo.csv", stringsAsFactors = F,
                 na.strings = c("NA", "", " ", "-"))

##### Prepare data #####
# Column file name
names(rawData)[names(rawData) %in% "QPCR01.06.2018.XLS.csv"] <- "fileName"

# Annoying spaces
rawData[,names(rawData)  %in% "fileName"] <- gsub(" ", "", rawData[,names(rawData) %in% "fileName"])
rawData[,names(rawData)  %in% "Target.SYBR"] <- gsub(" ", "", rawData[,names(rawData) %in% "Target.SYBR"])

# Remove inside headers
rawData <- rawData[rawData$Name != "Name",]

rawData[rawData$fileName %in% "QPCR03.05.2018_complete.csv", "fileName"] <- "QPCR03.05.2018.XLS.csv"  

rawData$fileName[rawData$fileName %in% "QPCR311-317.XLS.csv"] <- "QPCR09.05.2018.XLS.csv"

## Add melting curves infos when available
rawMeltData <- read.csv("./LorenzoRAW/MeltingCurves/MeltingCurvesLorenzo.csv",
                        na.strings = c("NA", "", " ", "-"))

names(rawMeltData)[names(rawMeltData) %in% "QPCR01.06.2018_MeltCurve.csv"] <- "fileName"

# Annoying spaces
rawMeltData[,names(rawMeltData)  %in% "fileName"] <- gsub(" ", "", rawMeltData[,names(rawMeltData) %in% "fileName"])

# Remove inside headers
rawMeltData <- rawMeltData[rawMeltData$Name != "Name",]

# Remove samples with no nbr of Tm value (extra lines)
rawMeltData <- rawMeltData[!is.na(rawMeltData$No..Tm.SYBR),]

# Correct fileName to match previous files
rawMeltData$fileName <- gsub("_MeltCurve", ".XLS", rawMeltData$fileName)

# Manual error : on plate 13.06.2018
rawMeltData$Name <- as.character(rawMeltData$Name)
rawMeltData[rawMeltData$fileName %in% "QPCR13.06.2018correct.XLS.csv" &
              rawMeltData$Pos %in% paste0("F", 7:12), "Name"] <- "CEWE_AA_0367"

rawData <- merge(rawData, 
                 rawMeltData[c("fileName", "Name", "Pos", "No..Tm.SYBR", "Tm.x..SYBR")],
                 by = c("fileName", "Name", "Pos"), all = T)

# Remove controls
rawData <- rawData[!rawData$Name %in% c("water", "NTC", "Noiseband","OFF", "Quantification", "SYBR","n/a" ) &
                     !rawData$Pos %in% c("Threshold level", "Baseline setting"),]

# Remove this plate as all samples were redone, better, later + no melt curve
rawData <- rawData[!rawData$fileName %in% "QPCR04.04.2018.XLS.csv",]

# No melt curves for them:
unique(rawData[is.na(rawData$No..Tm.SYBR),"fileName"])

# Remove stupid empty rows
rawData <- rawData[!is.na(rawData$fileName),]

# manual corrections on sample names
rawData$Name[rawData$Name %in% "359"] <- "ILWE_AA_0359"
rawData$Name[rawData$Name %in% "ILWE_AA_242"] <- "ILWE_AA_0242"
rawData$Name[rawData$Name %in% "CEWE_AA_0330#"] <- "CEWE_AA_0330"
rawData$Name[rawData$Name %in% "ILWE_AA_0,48"] <- "ILWE_AA_0348"
rawData$Name[rawData$Name %in% "ILWE_AA_0134"] <- "ILWE_AA_0379"
rawData[rawData$fileName %in% "QPCR07.06.2018part2.XLS.csv" &
          rawData$Pos %in% paste0("C", 4:12),"Name"] <- gsub("0410", "0409",
                                                             rawData[rawData$fileName %in% "QPCR07.06.2018part2.XLS.csv" &
                                                                       rawData$Pos %in% paste0("C", 4:12),"Name"])
rawData[rawData$Name %in% "CEWE_AA_0330#", "Name"] <- "CEWE_AA_0330"
rawData[rawData$Name %in% "CEWE_AA_0436*", "Name"] <- "CEWE_AA_0436"
 
rawData[rawData$Name %in% "CEWE_AA_0410" & rawData$Pos %in% c(paste0("D", 7:12)),"Name"] <- "ILWE_AA_0410"

# likely manual mistake
rawData[rawData$Pos %in% c("B4", "B5", "B6") &
          rawData$fileName %in% "QPCR14.06.2018.XLS.csv", "Name"] <- "CEWE_AA_0424"

# Add full name of sample (tissue + mouseID + eimeriaOrmouse primers + plate)
rawData$fullName <- paste0(rawData$Name, "_", rawData$Target.SYBR, "_",  rawData$fileName)

# Add tissue and Mouse_ID
x <- strsplit(as.character(rawData$Name), "_", 1)
rawData$tissue <- sapply( x, "[", 1)
rawData$Mouse_ID <- paste0("AA_", sapply( x, "[", 3))
rm(x)
##### END Prepare data #####

##### End cleaning ##### 
allSamples <- unique(data.frame(Mouse_ID = rawData$Mouse_ID,
                         tissue = rawData$tissue))

### 1. which ones are at least duplicates? remove the others
# triplicate = same file, same name, same target, same mean

# Remove samples with no delta Ct value for mice (means that only one sample worked)
failedDNA <- unique(data.frame(Mouse_ID = rawData[which(is.na(rawData$Ct.Mean.SYBR) & rawData$Target.SYBR == "mouse"), "Mouse_ID"],
                               tissue = rawData[which(is.na(rawData$Ct.Mean.SYBR) & rawData$Target.SYBR == "mouse"), "tissue"],
                               status = "failed"))

rawData <- rawData[-which(is.na(rawData$Ct.Mean.SYBR) & rawData$Target.SYBR == "mouse"),]

library(dplyr)

sumOKData <- rawData %>% 
  group_by(fullName) %>% 
  summarise(count = n()) %>% 
  data.frame()

OKData <- rawData[!rawData$fullName %in% 
                   sumOKData[sumOKData$count < 3, ]$fullName, ]

# 2. all samples with no CtMean from Eimeria (but with CtMean from mouse, cleaned before) are NEGATIVE
OKData <- OKData %>%
  group_by(fullName) %>%
  mutate(status = if_else(is.na(Ct.Mean.SYBR), "negative", "pending"))%>% 
  data.frame()

negative <- unique(data.frame(Mouse_ID = OKData[OKData$status == "negative", "Mouse_ID"],
                              tissue = OKData[OKData$status == "negative", "tissue"],
                              status = "negative"))

##### Calculate delta ct #####
calculateDeltaCt <- function(df){
  # Keep one value per plate per fullName (so per triplicate)
  df <- df[!duplicated(df[c("fullName")]),]
  ## Calculate deltaCt per plate
  sumDataMouse <- df[df$Target.SYBR %in% "mouse",]
  sumDataEimeria <- df[df$Target.SYBR %in% "eimeria",]
  mergedData <- merge(sumDataEimeria, sumDataMouse,
                      by = c("Mouse_ID", "fileName", "tissue", "Name"))
  mergedData <- unique(mergedData)
  mergedData$deltaCtMminusE <- as.numeric(as.character(mergedData$Ct.Mean.SYBR.y)) - 
    as.numeric(as.character(mergedData$Ct.Mean.SYBR.x)) # DeltaCt MOUSE minus EIMERIA
  return(mergedData)
}

OKData <- calculateDeltaCt(OKData)

hist(OKData$deltaCtMminusE, breaks = 100) # keep all above -6 :D

negative2 <- data.frame(Mouse_ID = OKData[OKData$deltaCtMminusE <= -6, "Mouse_ID"],
                        tissue = OKData[OKData$deltaCtMminusE <= -6, "tissue"],
                        status = "negative")
                        

positiveData <- OKData[OKData$deltaCtMminusE > -6,]

positiveSamples <- unique(data.frame(Mouse_ID = positiveData["Mouse_ID"],
                                    tissue = positiveData["tissue"],
                                    status = "positive"))

# Prevalence 2017
nrow(positiveSamples) / nrow(allSamples) *100

## Average technical replicates
finalData <- positiveData %>% 
  group_by(Name) %>%   
  summarise(count = n(), 
            deltaCt_MminusE = mean(deltaCtMminusE),
            sdDelta = sd(deltaCtMminusE)) %>% 
  data.frame()

# Final DF
x <- strsplit(as.character(finalData$Name), "_", 1)
finalData <- data.frame(Name = finalData$Name,
                        deltaCt_MminusE = finalData$deltaCt_MminusE,
                        tissue = sapply( x, "[", 1),
                        Mouse_ID = paste0("AA_", sapply( x, "[", 3)))
rm(x)

finalData$year <- 2017

finalDataClean <- finalData["Mouse_ID"]

# Add CEWE
finalDataClean <- merge(finalDataClean,
                        finalData[finalData$tissue %in% "CEWE", c("Mouse_ID", "deltaCt_MminusE")],
                           all.x = T)
names(finalDataClean)[names(finalDataClean) %in% "deltaCt_MminusE"] <- "delta_ct_cewe_MminusE"

# Add ILWE
finalDataClean <- merge(finalDataClean,
                           finalData[finalData$tissue %in% "ILWE", c("Mouse_ID", "deltaCt_MminusE")],
                           all.x = T)
names(finalDataClean)[names(finalDataClean) %in% "deltaCt_MminusE"] <- "delta_ct_ilwe_MminusE"

finalDataClean <- unique(finalDataClean)

# Add observer
finalDataClean$observer_qpcr <- "Lorenzo"

# Add year
finalDataClean$Year <- 2017
finalDataClean <- finalDataClean[-which(finalDataClean$Mouse_ID %in% "AA_NA"),]

# Write out
write.csv(finalDataClean, "../qPCR_2017.csv", row.names = F)
