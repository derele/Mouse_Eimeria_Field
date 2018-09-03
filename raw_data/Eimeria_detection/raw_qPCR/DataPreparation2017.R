# Back to empty
rm(rawData)
rm(rawMeltData)

# import data
rawData <- read.csv("./LorenzoRAW/CSVFiles/TotalLorenzo.csv", stringsAsFactors = F,
                 na.strings = c("NA", "", " ", "-"))

##### Clean data #####
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
rawData$Name[rawData$Name %in% "ILWE_AA_0,48"] <- "ILWE_AA_0348"
rawData$Name[rawData$Name %in% "ILWE_AA_0134"] <- "ILWE_AA_0379"
rawData[rawData$fileName %in% "QPCR07.06.2018part2.XLS.csv" &
          rawData$Pos %in% paste0("C", 4:12),"Name"] <- gsub("0410", "0409",
                                                             rawData[rawData$fileName %in% "QPCR07.06.2018part2.XLS.csv" &
                                                                       rawData$Pos %in% paste0("C", 4:12),"Name"])
rawData[rawData$Name %in% "CEWE_AA_0330#", "Name"] <- "CEWE_AA_0330"
rawData[rawData$Name %in% "CEWE_AA_0410" & rawData$Pos %in% c(paste0("D", 7:12)),"Name"] <- "ILWE_AA_0410"

# likely manual mistake
rawData[rawData$Pos %in% c("B4", "B5", "B6") &
          rawData$fileName %in% "QPCR14.06.2018.XLS.csv", "Name"] <- "CEWE_AA_0424"

# Remove samples with no delta Ct value for mice (means that only one sample worked)
rawData <- rawData[-which(is.na(rawData$Ct.Mean.SYBR) & rawData$Target.SYBR == "mouse"),]

# Add full name of sample (tissue + mouseID + eimeriaOrmouse primers + plate)
rawData$fullName <- paste0(rawData$Name, "_", rawData$Target.SYBR, "_",  rawData$fileName)

# Add tissue and Mouse_ID
x <- strsplit(as.character(rawData$Name), "_", 1)
rawData$tissue <- sapply( x, "[", 1)
rawData$Mouse_ID <- paste0("AA_", sapply( x, "[", 3))
rm(x)

##### End cleaning ##### 
allSamples <- unique(rawData$Name)

### 1. which ones are at least duplicates? remove the others
# triplicate = same file, same name, same target, same mean

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

##### Calculate delta ct #####
calculateDeltaCt <- function(df){
  # Keep one value per plate per fullName (so per triplicate)
  df <- df[!duplicated(df[c("fullName")]),]
  ## Calculate deltaCt per plate
  sumDataMouse <- df[df$Target.SYBR %in% "mouse",]
  sumDataEimeria <- df[df$Target.SYBR %in% "eimeria",]
  
  mergedData <- merge(sumDataEimeria, sumDataMouse, 
                      by = c("Mouse_ID", "fileName", "tissue"), all = T)
  
  mergedData$deltaCtMminusE <- NA
  
  mergedData$deltaCtMminusE[!is.na(as.numeric(mergedData$Ct.Mean.SYBR.y)) &
                              !is.na(as.numeric(mergedData$Ct.Mean.SYBR.x))] <- 
    as.numeric(mergedData$Ct.Mean.SYBR.y[!is.na(as.numeric(mergedData$Ct.Mean.SYBR.y)) &
                                           !is.na(as.numeric(mergedData$Ct.Mean.SYBR.x))]) - 
    as.numeric(mergedData$Ct.Mean.SYBR.x[!is.na(as.numeric(mergedData$Ct.Mean.SYBR.y)) &
                                           !is.na(as.numeric(mergedData$Ct.Mean.SYBR.x))])
  return(mergedData)
}

OKData <- calculateDeltaCt(OKData)

hist(OKData$deltaCtMminusE, breaks = 100) # keep all above -6 :D

OKData$status <- "negative"
OKData$status[OKData$deltaCtMminusE >= -6] <- "positive"

positiveData <- OKData[OKData$status %in% "positive",]

## Average technical replicates
finalData <- positiveData %>% 
  group_by(Name.x) %>%   
  summarise(count = n(), 
            deltaCt_MminusE = mean(deltaCtMminusE),
            sdDelta = sd(deltaCtMminusE)) %>% 
  data.frame()

positiveSample <- finalData$Name.x

# Prevalence 2017
length(positiveSample) / length(allSamples) *100

# Final DF
x <- strsplit(as.character(finalData$Name.x), "_", 1)
finalData <- data.frame(Name = finalData$Name.x,
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

# Write out
write.csv(finalDataClean, "../qPCR_2017.csv", row.names = F)

# ###########
# source("../../../R/functions/addPCRresults.R")
# source("../../../R/functions/addqPCRresults.R")
# source("../../../R/functions/addFlotationResults.R")
# 
# myFinal <- addPCRresults(finalData, pathtodata = "../Inventory_contents_all.csv")
# myFinal <- addFlotationResults(myFinal, pathtofinalOO = "../FINALOocysts2015to2017.csv",
#                                pathtolorenzodf = "../Eimeria_oocysts_2015&2017_Lorenzo.csv")$new
# 
# myFinal <- addqPCRresults(myFinal, pathtoqPCR2016 = "../qPCR_2016.csv", pathtoqPCR2017 = "../qPCR_2017.csv")
# 
# myFinal$year[is.na(myFinal$year)] <- myFinal$year.x[is.na(myFinal$year)]
# myFinal$year[is.na(myFinal$year)] <- myFinal$year.y[is.na(myFinal$year)]
# myFinal$year <- as.factor(myFinal$year)
# 
# summary(lm(OPG ~ delta_ct_cewe, myFinal))
# 
# ggplot(myFinal, aes(y = myFinal$delta_ct_ilwe, x = OPG+1)) +
#   scale_x_log10() +
#   geom_point(aes(col = year), size = 4) +
#   geom_smooth(method = "lm")
# 
# ggplot(myFinal, aes(y = myFinal$delta_ct_cewe, x = OPG+1)) +
#   scale_x_log10() +
#   geom_point(aes(col = year), size = 4) +
#   geom_smooth(method = "lm") 
# 
# # Plot 1 detection methods compared
# ggplot(myFinal, aes(x = PCRstatus, y = OPG + 1)) +
#   scale_y_log10() +
#   # geom_boxplot()+
#   geom_violin() +
#   geom_jitter(aes(col = year), width = .1, size = 2, alpha = .5) +
#   theme_bw()
# 
# # Plot 2 detection methods compared
# ggplot(myFinal, aes(x = qPCRstatus, y = OPG + 1)) +
#   scale_y_log10() +
#   geom_boxplot()+
#   geom_jitter(aes(col = year), width = .1, size = 2, alpha = .5) +
#   theme_bw()
# 
# # How to set up a limit of detection for qPCR
# ggplot(myFinal, aes(x = OPG > 0, y = delta_ct_cewe)) +
#   geom_boxplot()+
#   geom_point(aes(col = year), size = 2, alpha = .5) +
#   theme_bw()
# 
# ggplot(myFinal, aes(x = myFinal$PCRstatus, y = delta_ct_cewe)) +
#   geom_boxplot()+
#   geom_point(aes(col = year), size = 2, alpha = .5) +
#   theme_bw()
# 
# myFinal$FlotOrPcr <- "negative"
# myFinal$FlotOrPcr[myFinal$OPG > 0 | myFinal$PCRstatus %in% "positive"] <- "positive"
# 
# ggplot(myFinal, aes(x = myFinal$FlotOrPcr, y = delta_ct_cewe)) +
#   geom_violin()+
#   geom_point(aes(col = year), size = 2, alpha = .5) +
#   geom_hline(yintercept = 6) +
#   theme_bw()
# 
# ## Which samples have no qPCR?
# # 2017
# myFinal$Mouse_ID[!is.na(myFinal$delta_ct_cewe)]
# which(duplicated(myFinal$Mouse_ID))
# 





