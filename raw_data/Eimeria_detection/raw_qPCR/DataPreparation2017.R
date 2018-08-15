# import data
rawData <- read.csv("./LorenzoRAW/CSVFiles/TotalLorenzo.csv", stringsAsFactors = F,
                 na.strings = c("NA", "", " ", "-"))

##### Clean data #####
# Column file name
names(rawData)[names(rawData) == "QPCR01.06.2018.XLS.csv"] <- "fileName"

# Annoying spaces
rawData[,names(rawData)  == "fileName"] <- gsub(" ", "", rawData[,names(rawData) == "fileName"])
rawData[,names(rawData)  == "Target.SYBR"] <- gsub(" ", "", rawData[,names(rawData) == "Target.SYBR"])

# Remove inside headers
rawData <- rawData[rawData$Name != "Name",]

# Remove samples with no Ct value
rawData <- rawData[!is.na(rawData$Ct.SYBR),]

# Remove samples with no mean Ct (means that only one sample worked)
rawData <- rawData[!is.na(rawData$Ct.Mean.SYBR),]

# Remove controls (were used before)
rawData <- rawData[!rawData$Name %in% c("water", "NTC"),]

# manual correction
rawData$Name[rawData$Name == "359"] <- "ILWE_AA_0359"
rawData$Name[rawData$Name == "ILWE_AA_242"] <- "ILWE_AA_0242"
rawData$Name[rawData$Name == "ILWE_AA_0,48"] <- "ILWE_AA_0348"
rawData$Name[rawData$Name == "ILWE_AA_0134"] <- "ILWE_AA_0379"
rawData$fileName[rawData$fileName == "QPCR03.05.2018_complete.csv"] <- "QPCR03.05.2018.XLS.csv"  

# likely manual mistake
rawData[rawData$Pos %in% c("B4", "B5", "B6") & 
       rawData$fileName == "QPCR14.06.2018.XLS.csv", "Name"] <- "CEWE_AA_0424"

# Add full name of sample (tissue + mouseID + eimeriaOrmouse primers + plate)
rawData$fullName <- paste0(rawData$Name, "_", rawData$Target.SYBR, "_",  rawData$fileName)

# Add tissue and Mouse_ID
x <- strsplit(as.character(rawData$Name), "_", 1)
rawData$tissue <- sapply( x, "[", 1)
rawData$Mouse_ID <- paste0("AA_", sapply( x, "[", 3))
rm(x)

## Add melting curves infos when available
rawMeltData <- read.csv("./LorenzoRAW/MeltingCurves/MeltingCurvesLorenzo.csv",
                        na.strings = c("NA", "", " ", "-"))

names(rawMeltData)[names(rawMeltData) == "QPCR01.06.2018_MeltCurve.csv"] <- "fileName"

# Annoying spaces
rawMeltData[,names(rawMeltData)  == "fileName"] <- gsub(" ", "", rawMeltData[,names(rawMeltData) == "fileName"])

# Remove inside headers
rawMeltData <- rawMeltData[rawMeltData$Name != "Name",]

# Remove samples with no nbr of Tm value
rawMeltData <- rawMeltData[!is.na(rawMeltData$No..Tm.SYBR),]

# Correct fileName to match previous files
rawMeltData$fileName <- gsub("_MeltCurve", ".XLS", rawMeltData$fileName)

rawData <- merge(rawData, 
                 rawMeltData[c("fileName", "Name", "Pos", "No..Tm.SYBR")],
      by = c("fileName", "Name", "Pos"), all.x = T)

rm(rawMeltData)

##### End cleaning ##### 

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
  gg <- gg + geom_tile(aes(x = variable, y = tissue_target, fill = value), color="#7f7f7f")
  # custom fill colors
  gg <- gg + scale_fill_manual(values=c("grey", "green"))
  # squares
  gg <- gg + coord_equal()
  # no labels
  gg <- gg + labs(x=NULL, y=NULL)
  # remove some chart junk
  gg <- gg + theme_bw() + theme(panel.grid=element_blank(),
                                panel.border=element_blank(),
                                axis.text.x = element_text(angle = 45, hjust = 1, size=5) )
  return(table(dat_long$tissue_target, dat_long$value))
  gg
}

myTiles(rawData)

# Separate here failed and ok data
OKData <- rawData[!is.na(rawData$No..Tm.SYBR) &
                          rawData$No..Tm.SYBR != 0,]
myTiles(OKData)

# First,let's check OKData
## Which samples are complete (mouse+eimeria triplicate, low sd)

# CEWE_AA_0410_mouse_QPCR07.06.2018part2.XLS.csv 5 occurences, remove
OKData <- OKData[OKData$fullName != "CEWE_AA_0410_mouse_QPCR07.06.2018part2.XLS.csv",]
myTiles(OKData)

### 1. which ones are triplicates? remove the others
# triplicate = same file, same name, same target, same mean

library(dplyr)

sumOKData <- OKData %>% 
  group_by(fullName) %>% 
  summarise(count = n()) %>% 
  data.frame()

OKData <- OKData[!OKData$fullName %in% 
                   sumOKData[sumOKData$count != 3, ]$fullName, ]
myTiles(OKData)

# 2. how are the sd? keep < 3
OKData <- OKData[OKData$Ct.Dev..SYBR <= 3 ,]
myTiles(OKData)

# 3. remove duplicates on different plates (based on sd?)
OKData$nameTarget <- paste(OKData$Name, OKData$Target.SYBR)

sumOKData <- OKData %>% 
  group_by(nameTarget) %>% 
  summarise(isDup = length(nameTarget)) %>% 
  data.frame()

myTiles(OKData)

# 4. Eventually, mouse AND eimeria on the same plate otherwise remove
OKData$fileNameName <- paste(OKData$fileName, OKData$Name)

sumOKData <- OKData %>% 
  group_by(fileNameName) %>% 
  summarise(isBoth = length(unique(Target.SYBR))) %>% 
  data.frame()

OKData <- OKData[OKData$fileNameName %in% sumOKData$fileNameName[sumOKData$isBoth == 2],]

myTiles(OKData)

# No melt curves for:
# QPCR04.04.2018.XLS.csv but redonne twice later
# QPCR30.05.2018.XLS.csv --> check manually the look of the curves. Keep for later
# QPCR311-317.XLS.csv --> check manually the look of the curves. Keep for later

##### Calculate delta ct #####

calculateDeltqCt <- function(df){
  # Keep one value per plate per fullName (so per triplicate)
  df <- df[!duplicated(df[c("fullName")]),]
  
  ## Calculate deltaCt per plate
  sumDataMouse <- df[df$Target.SYBR == "mouse",]
  sumDataEimeria <- df[df$Target.SYBR == "eimeria",]
  
  mergedData <- merge(sumDataEimeria, sumDataMouse, 
                      by = c("Mouse_ID", "fileName", "tissue"))
  
  mergedData$deltaCt <- as.numeric(mergedData$Ct.Mean.SYBR.x) - as.numeric(mergedData$Ct.Mean.SYBR.y)
  return(mergedData)
}

finalData <- calculateDeltqCt(OKData)

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
                        finalData[finalData$tissue == "CEWE", c("Mouse_ID", "deltaCt")],
                           all.x = T)
names(finalDataClean)[names(finalDataClean) == "deltaCt"] <- "delta_ct_cewe"

# Add ILWE
finalDataClean <- merge(finalDataClean,
                           finalData[finalData$tissue == "ILWE", c("Mouse_ID", "deltaCt")],
                           all.x = T)
names(finalDataClean)[names(finalDataClean) == "deltaCt"] <- "delta_ct_ilwe"

# Add observer
finalDataClean$observer_qpcr <- "Lorenzo"

# Write out
write.csv(finalData, "../qPCR_2017.csv", row.names = F)

###########
source("../../../R/functions/addPCRresults.R")
source("../../../R/functions/addqPCRresults.R")
source("../../../R/functions/addFlotationResults.R")

myFinal <- addPCRresults(finalData, pathtodata = "../Inventory_contents_all.csv")
myFinal <- addFlotationResults(myFinal, pathtofinalOO = "../FINALOocysts2015to2017.csv",
                               pathtolorenzodf = "../Eimeria_oocysts_2015&2017_Lorenzo.csv")$new

myFinal <- addqPCRresults(myFinal, pathtoqPCR2016 = "../qPCR_2016.csv", pathtoqPCR2017 = "../qPCR_2017.csv")

myFinal$year[is.na(myFinal$year)] <- myFinal$year.x[is.na(myFinal$year)]
myFinal$year[is.na(myFinal$year)] <- myFinal$year.y[is.na(myFinal$year)]
myFinal$year <- as.factor(myFinal$year)

summary(lm(OPG ~ delta_ct_cewe, myFinal))

ggplot(myFinal, aes(y = myFinal$delta_ct_ilwe, x = OPG)) +
  geom_point(aes(col = year), size = 4) +
  geom_smooth(method = "lm")

# Plot 1 detection methods compared
ggplot(myFinal, aes(x = PCRstatus, y = OPG + 1)) +
  scale_y_log10() +
  # geom_boxplot()+
  geom_violin() +
  geom_jitter(aes(col = year), width = .1, size = 2, alpha = .5) +
  theme_bw()

# Plot 2 detection methods compared
ggplot(myFinal, aes(x = qPCRstatus, y = OPG + 1)) +
  scale_y_log10() +
  geom_boxplot()+
  geom_jitter(aes(col = year), width = .1, size = 2, alpha = .5) +
  theme_bw()

# How to set up a limit of detection for qPCR
ggplot(myFinal, aes(x = OPG > 0, y = delta_ct_cewe)) +
  geom_boxplot()+
  geom_point(aes(col = year), size = 2, alpha = .5) +
  theme_bw()

ggplot(myFinal, aes(x = myFinal$PCRstatus, y = delta_ct_cewe)) +
  geom_boxplot()+
  geom_point(aes(col = year), size = 2, alpha = .5) +
  theme_bw()

myFinal$FlotOrPcr <- "negative"
myFinal$FlotOrPcr[myFinal$OPG > 0 | myFinal$PCRstatus == "positive"] <- "positive"

ggplot(myFinal, aes(x = myFinal$FlotOrPcr, y = delta_ct_cewe)) +
  geom_violin()+
  geom_point(aes(col = year), size = 2, alpha = .5) +
  geom_hline(yintercept = 6) +
  theme_bw()

## Which samples have no qPCR?
# 2017
myFinal$Mouse_ID[!is.na(myFinal$delta_ct_cewe)]
which(duplicated(myFinal$Mouse_ID))

