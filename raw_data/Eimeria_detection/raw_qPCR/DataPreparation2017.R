rawMeltData <- read.csv("./LorenzoRAW/CSVFiles/TotalLorenzo.csv", stringsAsFactors = F,
                 na.strings = c("NA", "", " ", "-"))

##### Clean data #####
# Column file name
names(rawMeltData)[names(rawMeltData) == "QPCR01.06.2018.XLS.csv"] <- "fileName"

# Annoying spaces
rawMeltData[,names(rawMeltData)  == "fileName"] <- gsub(" ", "", rawMeltData[,names(rawMeltData) == "fileName"])
rawMeltData[,names(rawMeltData)  == "Target.SYBR"] <- gsub(" ", "", rawMeltData[,names(rawMeltData) == "Target.SYBR"])

# Remove inside headers
rawMeltData <- rawMeltData[rawMeltData$Name != "Name",]

# Remove samples with no Ct value
rawMeltData <- rawMeltData[!is.na(rawMeltData$Ct.SYBR),]

# Remove samples with no mean Ct (means that only one sample worked)
rawMeltData <- rawMeltData[!is.na(rawMeltData$Ct.Mean.SYBR),]

# Remove controls (were used before)
rawMeltData <- rawMeltData[!rawMeltData$Name %in% c("water", "NTC"),]

# manual correction
rawMeltData$Name[rawMeltData$Name == "359"] <- "ILWE_AA_0359"
rawMeltData$Name[rawMeltData$Name == "ILWE_AA_242"] <- "ILWE_AA_0242"
rawMeltData$Name[rawMeltData$Name == "ILWE_AA_0,48"] <- "ILWE_AA_0348"
rawMeltData$Name[rawMeltData$Name == "ILWE_AA_0134"] <- "ILWE_AA_0379"
# likely manual mistake
rawMeltData[rawMeltData$Pos %in% c("B4", "B5", "B6") & 
       rawMeltData$fileName == "QPCR14.06.2018.XLS.csv", "Name"] <- "CEWE_AA_0424"

# Add full name of sample (tissue + mouseID + eimeriaOrmouse primers + plate)
rawMeltData$fullName <- paste0(rawMeltData$Name, "_", rawMeltData$Target.SYBR, "_",  rawMeltData$fileName)

# Add tissue and Mouse_ID
x <- strsplit(as.character(rawMeltData$Name), "_", 1)

rawMeltData$tissue <- sapply( x, "[", 1)

rawMeltData$Mouse_ID <- paste0("AA_", sapply( x, "[", 3))

##### Select correct data #####

# How many values per plate per fullname? (should be 3 max)
library(dplyr)
keepOnlyTriplicates <- rawMeltData %>%
  group_by(fullName)%>%
  count()

paste0("Problem with ", pull(keepOnlyTriplicates[keepOnlyTriplicates$n > 3,c("fullName")]))

rawMeltData <- rawMeltData[!rawMeltData$fullName %in%
                     pull(
                       keepOnlyTriplicates[keepOnlyTriplicates$n > 3,
                                           c("fullName")]),]

# Keep one value per plate per fullName (so per triplicate)
sumData <- rawMeltData[!duplicated(rawMeltData[c("fullName")]),]

# Choose what to do when several plates (choose based on sd or date?)
# in case of repeated samples, take the first one that worked (sd < 3)
# "Checking the files (not all), Lorenzo made new attempts when the variability 
# among replicates was too high (Sd >3) maybe would be nice to include the Sd 
# value for each sample as a criteria for selection of the data... 
# He replicate complete plates, so in most of the cases he include samples even 
# when in the previous attempt they had a good result"

sumData$IdTargetTissue <- paste(sumData$Mouse_ID, sumData$Target.SYBR, sumData$tissue)

duplicatedData <- unique(sumData$IdTargetTissue[duplicated(sumData$IdTargetTissue)])
                    
# Split in 2 to work only on duplicated DF     
myData1 <- sumData[!sumData$IdTargetTissue %in% duplicatedData, ]
myData2 <- sumData[sumData$IdTargetTissue %in% duplicatedData, ]

# keep the lower sd for each sample "IdTargetTissue"
myData3 <- data.frame()
for (i in unique(myData2$IdTargetTissue)){
  sub <- myData2[myData2$IdTargetTissue == i,]
  myData3 <- rbind(myData3, sub[sub$Ct.Dev..SYBR %in% min(sub$Ct.Dev..SYBR),])
}

# Unique values here:  
sumData <- rbind(myData1, myData3)

## Calculate deltaCt per plate

sumDataMouse <- sumData[sumData$Target.SYBR == "mouse",]
sumDataEimeria <- sumData[sumData$Target.SYBR == "eimeria",]

mergedData <- merge(sumDataEimeria, sumDataMouse, by = c("Mouse_ID", "fileName", "tissue"))

mergedData$deltaCt <- as.numeric(mergedData$Ct.Mean.SYBR.x) - as.numeric(mergedData$Ct.Mean.SYBR.y)

library(ggplot2)
ggplot(mergedData, aes(x = mergedData$Mouse_ID, y = mergedData$deltaCt)) +
  geom_point(aes(col = mergedData$tissue))

###### Calculation of LOD (mean + 2 standard deviations of the negative controls) #####

rawMeltData <- read.csv("./LorenzoRAW/MeltingCurves/MeltingCurvesLorenzo.csv", stringsAsFactors = F,
                    na.strings = c("NA", "", " ", "-"))

##### Clean data #####
# Column file name
names(rawMeltData)[names(rawMeltData) == "QPCRmc06.04.2018.XLS.csv"] <- "fileName"

# Annoying spaces
rawMeltData[,names(rawMeltData)  == "fileName"] <- gsub(" ", "", rawMeltData[,names(rawMeltData) == "fileName"])
rawMeltData[,names(rawMeltData)  == "Target.SYBR"] <- gsub(" ", "", rawMeltData[,names(rawMeltData) == "Target.SYBR"])

# Remove inside headers
rawMeltData <- rawMeltData[rawMeltData$Name != "Name",]

# Remove samples with no nbr of Tm value
rawMeltData <- rawMeltData[!is.na(rawMeltData$No..Tm.SYBR),]

# Correct fileName to match previous files
rawMeltData$fileName <- gsub("QPCRmc", "QPCR", rawMeltData$fileName)




# Remove samples with no mean Ct (means that only one sample worked)
rawMeltData <- rawMeltData[!is.na(rawMeltData$Ct.Mean.SYBR),]

# KEEP controls
ControlData <- rawMeltData[rawMeltData$Name %in% c("water", "NTC"),]

# Add full name of sample (control + eimeriaOrmouse primers + plate)
ControlData$fullName <- paste0(ControlData$Name, "_", ControlData$Target.SYBR, "_",  ControlData$fileName)

##### Select correct data #####

# How many values per plate per fullname? (should be 3 max)
library(dplyr)
keepOnlyTriplicates <- ControlData %>%
  group_by(fullName)%>%
  count()

# Keep one value per plate per fullName (so per triplicate)
sumData <- ControlData[!duplicated(ControlData[c("fullName")]),]

## Calculate deltaCt per plate per control type

sumDataMouse <- sumData[sumData$Target.SYBR == "mouse",]
sumDataEimeria <- sumData[sumData$Target.SYBR == "eimeria",]

mergedData <- merge(sumDataEimeria, sumDataMouse, by = c("fileName", "Name"))

mergedData$deltaCt <- as.numeric(mergedData$Ct.Mean.SYBR.x) - as.numeric(mergedData$Ct.Mean.SYBR.y)

# Calculation of LOD (mean + 2 standard deviations of the negative controls)

mean(mergedData$deltaCt) + 2 * sd(mergedData$deltaCt)

## NOT ENOUGH CONTROLS!!!

