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

# Remove samples with no nbr of Tm value
rawMeltData <- rawMeltData[!is.na(rawMeltData$No..Tm.SYBR),]

# Correct fileName to match previous files
rawMeltData$fileName <- gsub("_MeltCurve", ".XLS", rawMeltData$fileName)

# Manual error : on plate 13.06.2018
rawMeltData$Name <- as.character(rawMeltData$Name)
rawMeltData[rawMeltData$fileName %in% "QPCR13.06.2018correct.XLS.csv" &
              rawMeltData$Pos %in% paste0("F", 7:12), "Name"] <- "CEWE_AA_0367"

rawData <- merge(rawData, 
                 rawMeltData[c("fileName", "Name", "Pos", "No..Tm.SYBR")],
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
rawData[rawData$fileName %in% "QPCR07.06.2018part2.XLS.csv" &
          rawData$Pos %in% paste0("D", 7:12),"Name"] <- "ILWE_AA_0410" 
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

## heatmap to follow
library(ggplot2)

myTiles <- function(df){
  dat_long <- data.frame(tissue_target = paste(rawData$tissue, rawData$Target.SYBR), 
                         variable = rawData$Mouse_ID,
                         value = 0)
  
  dat_long <- unique(dat_long)
  dat_long[paste(dat_long$tissue_target, dat_long$variable) %in% 
             paste(df$tissue, df$Target.SYBR, df$Mouse_ID), "value"] <- 1
  
  # discrete vs continuous
  dat_long$value <- factor(dat_long$value)
  
  gg <- ggplot(dat_long)
  # fill + legend, gray border
  gg <- gg + geom_tile(aes(x = variable, y = tissue_target, fill = value),
                       color="#7f7f7f")
  # custom fill colors
  gg <- gg + scale_fill_manual(values=c("grey", "green"))
  # no labels
  gg <- gg + labs(x=NULL, y=NULL)
  # remove some chart junk
  gg <- gg + theme_bw() + theme(panel.grid=element_blank(),
                                panel.border=element_blank(),
                                axis.text.x = element_text(angle = 45, hjust = 1, size=3) )
  return(list(gg, table(dat_long$tissue_target, dat_long$value)))
}

myTiles(rawData)

# Separate here failed and ok data
# OKData <- rawData[!is.na(rawData$No..Tm.SYBR) &
#                           rawData$No..Tm.SYBR != 0,]
# myTiles(OKData)

##### Clean data #####
## 1. Which samples are complete (mouse+eimeria triplicate, low sd)
# triplicate = same file, same name, same target, same mean

library(dplyr)

sumOKData <- rawData %>% 
  group_by(fullName) %>% 
  summarise(count = n()) %>% 
  data.frame()

OKData <- rawData[!rawData$fullName %in% 
                   sumOKData[sumOKData$count < 3 , ]$fullName, ]
myTiles(OKData)

# 2. how are the sd? keep < 3
badSd <- OKData[!is.na(OKData$Ct.Dev..SYBR) & OKData$Ct.Dev..SYBR > 3 ,]
OKData <- OKData[OKData$Ct.Dev..SYBR <= 3 ,]
myTiles(OKData)

# 3. Keep plates where at least Mouse primers worked (our amplification control)
mouseData <- OKData[OKData$Target.SYBR %in% "mouse" & OKData$No..Tm.SYBR %in% 1,]
eimeriaData <- OKData[OKData$Target.SYBR %in% "eimeria",]

OKData <- rbind(mouseData, eimeriaData)

# 4. Calculate delta ct

calculateDeltaCt <- function(df){
  # Keep one value per plate per fullName (so per triplicate)
  df <- df[!duplicated(df[c("fullName")]),]
  
  ## Calculate deltaCt per plate
  sumDataMouse <- df[df$Target.SYBR %in% "mouse",]
  sumDataEimeria <- df[df$Target.SYBR %in% "eimeria",]
  
  mergedData <- merge(sumDataEimeria, sumDataMouse, 
                      by = c("Mouse_ID", "fileName", "tissue", "Name"))
  
  mergedData$deltaCt <- as.numeric(mergedData$Ct.Mean.SYBR.x) - as.numeric(mergedData$Ct.Mean.SYBR.y)
  return(mergedData)
}

OKData <- calculateDeltaCt(OKData)

# 5. Eventually, remove duplicates on different plates (based on sd)
sumOKData <- OKData %>% 
  group_by(Name) %>% 
  summarise(isDup = length(Name)) %>% 
  data.frame()

haveDuplicates <- OKData[OKData$Name %in% sumOKData$Name[sumOKData$isDup > 1],]

# Test: when several replicates, are they all going in the same direction (deltact > 6 or <6 for all)? 

testContradcition <- haveDuplicates %>% 
  group_by(Name) %>% 
  mutate(isContradictorybetweenDuplicates = length(table(deltaCt > 6))) %>% 
  data.frame()
table(testContradcition$isContradictorybetweenDuplicates) 
# hum, 29 tries contradictory. Check manually.

ggplot(haveDuplicates, aes(No..Tm.SYBR.x, deltaCt)) +
  geom_boxplot() +
  geom_point() +
  theme_bw() +
  scale_y_continuous(breaks = 0:20)
haveDuplicates$Mouse_ID


#### -> here i am lost
testContradcition$isContradictorybetweenDuplicates # perfect, always in the same sense :)

# we want to keep the replicate that minimize both sd

test <- toCheck[toCheck$Name == toCheck$Name[1],]

# per group
keepBestPerSample <- function(x){
  X = x[x$Ct.Dev..SYBR.x == min(x$Ct.Dev..SYBR.x),] # mouse minimal sd
  Y = x[x$Ct.Dev..SYBR.y == min(x$Ct.Dev..SYBR.y),] # eimeria minimal sd
  if(X$fileName == Y$fileName){
    unique <- rbind(unique, X)
  } else {
    unique <- rbind(unique, X, Y)
  }
  return(unique)
}

i = unique(toCheck$Name)[1]
x = toCheck[toCheck$Name %in% i,]

keepBestPerSample(x)

unique <- data.frame()

for (i in unique(toCheck$Name)){
  x = toCheck[toCheck$Name %in% i,]
  print(x)
  # unique = rbind(unique, keepBestPerSample(x))
}





by(data = toCheck, INDICES = toCheck$Name, function(x){keepBestPerSample(x)})

onEssaye <- toCheck %>% 
  group_by(Name) %>% 
  mutate(isContradictorybetweenDuplicates = length(table(test$deltaCt > 6))) %>% 
  data.frame()
lapply(min, toCheck$Ct.Dev..SYBR.x, toCheck$Ct.Dev..SYBR.y, min)

myTiles(OKData)

# Is Tm in Eimeria (=there should be DNA) linked to deltaCt...?
ggplot(finalData, aes(x = finalData$No..Tm.SYBR.x, y = deltaCt)) +
  geom_boxplot() +
  geom_jitter() +
  theme_bw()

library(ggplot2)
ggplot(finalData, aes(x = finalData$tissue, y = finalData$deltaCt)) +
  geom_violin() +
  geom_jitter(aes(col = finalData$tissue), size = 3) +
  theme_bw()

ggplot(finalData, aes(x = finalData$deltaCt)) +
  geom_histogram(aes(y=..density..), bins = 40) + 
  geom_density(aes(y=..density..)) +
  theme_bw()

finalData$year <- 2017

finalDataClean <- finalData["Mouse_ID"]
# Add CEWE
finalDataClean <- merge(finalDataClean,
                        finalData[finalData$tissue %in% "CEWE", c("Mouse_ID", "deltaCt")],
                           all.x = T)
names(finalDataClean)[names(finalDataClean) %in% "deltaCt"] <- "delta_ct_cewe"

# Add ILWE
finalDataClean <- merge(finalDataClean,
                           finalData[finalData$tissue %in% "ILWE", c("Mouse_ID", "deltaCt")],
                           all.x = T)
names(finalDataClean)[names(finalDataClean) %in% "deltaCt"] <- "delta_ct_ilwe"

# Add observer
finalDataClean$observer_qpcr <- "Lorenzo"

# Write out
write.csv(finalData, "../qPCR_2017.csv", row.names = F)

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

