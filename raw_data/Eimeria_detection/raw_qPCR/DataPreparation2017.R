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
# likely manual mistake
rawData[rawData$Pos %in% c("B4", "B5", "B6") & 
       rawData$fileName == "QPCR14.06.2018.XLS.csv", "Name"] <- "CEWE_AA_0424"

# Add full name of sample (tissue + mouseID + eimeriaOrmouse primers + plate)
rawData$fullName <- paste0(rawData$Name, "_", rawData$Target.SYBR, "_",  rawData$fileName)

# Add tissue and Mouse_ID
x <- strsplit(as.character(rawData$Name), "_", 1)

rawData$tissue <- sapply( x, "[", 1)

rawData$Mouse_ID <- paste0("AA_", sapply( x, "[", 3))

## Add melting curves infos

rawMeltData <- read.csv("./LorenzoRAW/MeltingCurves/MeltingCurvesLorenzo.csv", stringsAsFactors = F,
                        na.strings = c("NA", "", " ", "-"))

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

rawData <- merge(rawData, rawMeltData[c("fileName", "Name", "Pos", "No..Tm.SYBR")],
      by = c("fileName", "Name", "Pos"))

# Separate here positive and negative data
negativeData <- rawData[rawData$No..Tm.SYBR == "0",]
positiveData <- rawData[rawData$No..Tm.SYBR == "1",]

##### Select correct data and calculate delta ct #####

selectQPCRfun <- function(df){
  # How many values per plate per fullname? (should be 3 max)
  library(dplyr)
  keepOnlyTriplicates <- df %>%
    group_by(fullName)%>%
    count()
  
  paste0("Problem with ", pull(keepOnlyTriplicates[keepOnlyTriplicates$n > 3,c("fullName")]))
  
  df <- df[!df$fullName %in%
             pull(
               keepOnlyTriplicates[keepOnlyTriplicates$n > 3,
                                   c("fullName")]),]
  
  # Keep one value per plate per fullName (so per triplicate)
  sumData <- df[!duplicated(df[c("fullName")]),]
  
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
  
  return(mergedData)
}

finalPos <- selectQPCRfun(positiveData)
finalPos$isPos <- "positive"
finalNeg <- selectQPCRfun(negativeData)
finalNeg$isPos <- "negative"

finalDF <- rbind(finalPos, finalNeg)
finalDF$deltaCt

library(ggplot2)
ggplot(finalDF, aes(x = finalDF$Mouse_ID, y = finalDF$deltaCt)) +
  geom_point(aes(col = finalDF$isPos, pch = finalDF$tissue), size = 3) +
  theme_bw()

###### Calculation of LOD (mean + 2 standard deviations of the negative controls) #####
mean(finalNeg$deltaCt) + 2 * sd(finalNeg$deltaCt)

write.csv(finalPos, "../qPCR_2017.csv", row.names = F)
