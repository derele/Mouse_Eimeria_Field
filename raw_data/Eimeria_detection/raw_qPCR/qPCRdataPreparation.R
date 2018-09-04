### Prepare and homogenize qPCR data.
### Alice Balard
### So far, 2016 and 2017 samples

## Function to calculate deltaCtMminusE
calculateDeltaCt <- function(df, mergeBy){
  sumDataMouse <- df[df$Target.SYBR %in% "mouse",]
  sumDataEimeria <- df[df$Target.SYBR %in% "eimeria",]
  mergedData <- merge(sumDataEimeria, sumDataMouse,
                      by = mergeBy)
  mergedData <- unique(mergedData)
  mergedData$deltaCt_MminusE <- as.numeric(as.character(mergedData$Ct.Mean.SYBR.y)) - 
    as.numeric(as.character(mergedData$Ct.Mean.SYBR.x)) # DeltaCt MOUSE minus EIMERIA
  return(mergedData)
}

## Function to homogenize the final Df
myqPCRformat <- function(df){
  # which part is infected
  df$qPCRsummary[df$delta_ct_cewe_MminusE <= -6 & df$delta_ct_ilwe_MminusE <= -6] <- "non infected"
  df$qPCRsummary[df$delta_ct_cewe_MminusE > -6 & df$delta_ct_ilwe_MminusE <= -6] <- "infected cecum"
  df$qPCRsummary[df$delta_ct_cewe_MminusE <= -6 & df$delta_ct_ilwe_MminusE > -6] <- "infected ileum"
  df$qPCRsummary[df$delta_ct_cewe_MminusE > -6 & df$delta_ct_ilwe_MminusE > -6 & 
                   df$delta_ct_cewe_MminusE > df$delta_ct_ilwe_MminusE] <- "cecum stronger"
  df$qPCRsummary[df$delta_ct_cewe_MminusE > -6 & df$delta_ct_ilwe_MminusE > -6 & 
                   df$delta_ct_cewe_MminusE < df$delta_ct_ilwe_MminusE] <- "ileum stronger"
  df$qPCRsummary[df$delta_ct_cewe_MminusE > -6 & is.na(df$delta_ct_ilwe_MminusE)] <- "infected cecum [ileum.NA]"
  df$qPCRsummary[is.na(df$delta_ct_cewe_MminusE) & df$delta_ct_ilwe_MminusE > -6] <- "infected ileum [cecum.NA]"
  df$qPCRsummary[df$delta_ct_cewe_MminusE <= -6 & is.na(df$delta_ct_ilwe_MminusE)] <- "non infected cecum [ileum.NA]"
  df$qPCRsummary[is.na(df$delta_ct_cewe_MminusE) & df$delta_ct_ilwe_MminusE <= -6] <- "non infected ileum [cecum.NA]"
  
  # Infected or not? general status
  df$qPCRstatus <- "positive"
  df$qPCRstatus[is.na(df$qPCRsummary)] <- NA
  df$qPCRstatus[df$qPCRsummary %in% "non infected"] <- "negative"
  
  # and keep the infected segment value OR the higher value OR the only value (in some cases both were not tested)
  df$delta_ct_MminusE[df$qPCRsummary %in% c("infected cecum", "cecum stronger", "infected cecum [ileum.NA]")] <- 
    df$delta_ct_cewe_MminusE[df$qPCRsummary %in% c("infected cecum", "cecum stronger", "infected cecum [ileum.NA]")] 
  df$delta_ct_MminusE[df$qPCRsummary %in% c("infected ileum", "ileum stronger", "infected ileum [cecum.NA]")] <- 
    df$delta_ct_ilwe_MminusE[df$qPCRsummary %in% c("infected ileum", "ileum stronger", "infected ileum [cecum.NA]")] 
  # Set floor values for negative samples
  df$delta_ct_MminusE[is.na(df$delta_ct_MminusE)] <- -6
  # output always the same format
  data.frame(Mouse_ID = df$Mouse_ID,
             delta_ct_MminusE = df$delta_ct_MminusE,
             delta_ct_ilwe_MminusE = df$delta_ct_ilwe_MminusE,
             delta_ct_cewe_MminusE = df$delta_ct_cewe_MminusE,
             qPCRsummary = df$qPCRsummary,
             qPCRstatus = df$qPCRstatus,
             observer = df$observer_qpcr,
             year = df$year)
}

###### 2016 ######
qpcr2016 <- read.csv("MertRaw_2016/qPCR_2016.csv")
qpcr2016$year <- 2016
names(qpcr2016)[names(qpcr2016) %in% "sample.ID"] <- "Mouse_ID"

###### Enas data
qpcr2016_Enas <- qpcr2016[qpcr2016$observer_qpcr %in% "Enas",]  

# These are supposed to be positive samples (based on PCR)
# I assumed in this case it is deltact = CtMouse -CtEimeria
names(qpcr2016_Enas)[
  names(qpcr2016_Enas) %in% c("delta_ct_cewe", "delta_ct_ilwe")] <- 
  c("delta_ct_cewe_MminusE", "delta_ct_ilwe_MminusE")

# output final
qpcr2016_final_Enas <- myqPCRformat(qpcr2016_Enas)

###### Mert data
qpcr2016_Mert <- qpcr2016[qpcr2016$observer_qpcr %in% "Mert",]  
# Here deltaCT = ct eimeria - ct mouse. If high infection, low deltaCT
# Let's turn it around for higher value = higher load (so deltaCT = ct mouse - ct eimeria)
qpcr2016_Mert$delta_ct_cewe_MminusE <- - qpcr2016_Mert$delta_ct_cewe
qpcr2016_Mert$delta_ct_ilwe_MminusE <- - qpcr2016_Mert$delta_ct_ilwe

# output final
qpcr2016_final_Mert <- myqPCRformat(qpcr2016_Mert)

## Add Victor's values!
qpcr2016_Victor <- read.csv("VictorRaw_2016/Eimeria_qPCR_Tissue_220818.csv")
qpcr2016_MELT_Victor <- read.csv("VictorRaw_2016/Eimeria_qPCR_Tissue_220818_Melting.csv")

# Correct name
qpcr2016_Victor[qpcr2016_Victor$Name %in% "ILWE_AA_140", "Name"] <- "ILWE_AA_0140"

# Merge by Name and Pos
qpcr2016_Victor <- merge(qpcr2016_Victor, qpcr2016_MELT_Victor, by = c("Name", "Pos"))

# Remove useless lines
qpcr2016_Victor <- qpcr2016_Victor[!is.na(qpcr2016_Victor$Ct.Mean.SYBR) & !qpcr2016_Victor$Name %in% "NTC_Mouse", ]

# Add tissue and Mouse_ID
x <- strsplit(as.character(qpcr2016_Victor$Name), "_", 1)
qpcr2016_Victor$tissue <- sapply( x, "[", 1)
qpcr2016_Victor$Mouse_ID <- paste0("AA_", sapply( x, "[", 3))
rm(x)

# Calculate deltact
qpcr2016_Victor <- calculateDeltaCt(qpcr2016_Victor, 
                                   mergeBy = c("Mouse_ID", "tissue", "Name"))

qpcr2016_Victor <- unique(qpcr2016_Victor[c("Mouse_ID", "tissue", "deltaCt_MminusE")])

# Pre-formate
tempCEWE <- qpcr2016_Victor[qpcr2016_Victor$tissue == "CEWE",]
names(tempCEWE)[names(tempCEWE) %in% "deltaCt_MminusE"] <- "delta_ct_cewe_MminusE"

tempILWE <- qpcr2016_Victor[qpcr2016_Victor$tissue == "ILWE",]
names(tempILWE)[names(tempILWE) %in% "deltaCt_MminusE"] <- "delta_ct_ilwe_MminusE"

qpcr2016_Victor <- merge(tempCEWE, tempILWE, by = c("Mouse_ID"), all = T)

qpcr2016_Victor$year <- 2016
qpcr2016_Victor$observer_qpcr <- "Victor"

# output final
qpcr2016_final_Victor <- myqPCRformat(qpcr2016_Victor)

## Add Tabea's values!
qpcr2016_Tabea <- read.csv("TabeaRaw_2016/Other_Rodents_120718_3.csv")
qpcr2016_MELT_Tabea <- read.csv("TabeaRaw_2016/Other_Rodents_120718_Melting_3.csv")

# keep mus only
qpcr2016_Tabea <- qpcr2016_Tabea[grep("AA_", qpcr2016_Tabea$Name),]
qpcr2016_MELT_Tabea <- qpcr2016_MELT_Tabea[grep("AA_", qpcr2016_MELT_Tabea$Name),]

# Merge by Name and Pos
qpcr2016_Tabea <- merge(qpcr2016_Tabea, qpcr2016_MELT_Tabea, by = c("Name", "Pos"))

# replace with correct naming CeWe and Ilwe
qpcr2016_Tabea$Name <- gsub("CeWe ","CEWE_", qpcr2016_Tabea$Name)
qpcr2016_Tabea$Name <- gsub("ILWe ","ILWE_", qpcr2016_Tabea$Name)
qpcr2016_Tabea$Name <- gsub("IlWe ","ILWE_", qpcr2016_Tabea$Name)

# Add tissue and Mouse_ID
x <- strsplit(as.character(qpcr2016_Tabea$Name), "_", 1)
qpcr2016_Tabea$tissue <- sapply( x, "[", 1)
qpcr2016_Tabea$Mouse_ID <- paste0("AA_", sapply( x, "[", 3))
rm(x)

# Calculate deltact
qpcr2016_Tabea <- calculateDeltaCt(qpcr2016_Tabea, 
                                    mergeBy = c("Mouse_ID", "tissue", "Name"))

qpcr2016_Tabea <- unique(qpcr2016_Tabea[c("Mouse_ID", "tissue", "deltaCt_MminusE")])

# Pre-formate
tempCEWE <- qpcr2016_Tabea[qpcr2016_Tabea$tissue == "CEWE",]
names(tempCEWE)[names(tempCEWE) %in% "deltaCt_MminusE"] <- "delta_ct_cewe_MminusE"

tempILWE <- qpcr2016_Tabea[qpcr2016_Tabea$tissue == "ILWE",]
names(tempILWE)[names(tempILWE) %in% "deltaCt_MminusE"] <- "delta_ct_ilwe_MminusE"

qpcr2016_Tabea <- merge(tempCEWE, tempILWE, by = c("Mouse_ID"), all = T)

qpcr2016_Tabea$year <- 2016
qpcr2016_Tabea$observer_qpcr <- "Tabea"

# output final
qpcr2016_final_Tabea <- myqPCRformat(qpcr2016_Tabea)

###### 2017 ######

###### Lorenzo data

# import data
qpcr2017_Lorenzo <- read.csv("LorenzoRAW_2017/CSVFiles/TotalLorenzo.csv", stringsAsFactors = F,
                             na.strings = c("NA", "", " ", "-"))

##### Prepare data #####
# Column file name
names(qpcr2017_Lorenzo)[names(qpcr2017_Lorenzo) %in% "QPCR01.06.2018.XLS.csv"] <- "fileName"

# Annoying spaces
qpcr2017_Lorenzo[,names(qpcr2017_Lorenzo)  %in% "fileName"] <- 
  gsub(" ", "", qpcr2017_Lorenzo[,names(qpcr2017_Lorenzo) %in% "fileName"])
qpcr2017_Lorenzo[,names(qpcr2017_Lorenzo)  %in% "Target.SYBR"] <- 
  gsub(" ", "", qpcr2017_Lorenzo[,names(qpcr2017_Lorenzo) %in% "Target.SYBR"])

# Remove inside headers
qpcr2017_Lorenzo <- qpcr2017_Lorenzo[qpcr2017_Lorenzo$Name != "Name",]

qpcr2017_Lorenzo[qpcr2017_Lorenzo$fileName %in% "QPCR03.05.2018_complete.csv", "fileName"] <- "QPCR03.05.2018.XLS.csv"  

qpcr2017_Lorenzo$fileName[qpcr2017_Lorenzo$fileName %in% "QPCR311-317.XLS.csv"] <- "QPCR09.05.2018.XLS.csv"

## Add melting curves infos when available
qpcr2017_MELT_Lorenzo <- read.csv("LorenzoRAW_2017/MeltingCurves/MeltingCurvesLorenzo.csv",
                        na.strings = c("NA", "", " ", "-"))

names(qpcr2017_MELT_Lorenzo)[names(qpcr2017_MELT_Lorenzo) %in% "QPCR01.06.2018_MeltCurve.csv"] <- "fileName"

# Annoying spaces
qpcr2017_MELT_Lorenzo[,names(qpcr2017_MELT_Lorenzo)  %in% "fileName"] <- gsub(" ", "", qpcr2017_MELT_Lorenzo[,names(qpcr2017_MELT_Lorenzo) %in% "fileName"])

# Remove inside headers
qpcr2017_MELT_Lorenzo <- qpcr2017_MELT_Lorenzo[qpcr2017_MELT_Lorenzo$Name != "Name",]

# Remove samples with no nbr of Tm value (extra lines)
qpcr2017_MELT_Lorenzo <- qpcr2017_MELT_Lorenzo[!is.na(qpcr2017_MELT_Lorenzo$No..Tm.SYBR),]

# Correct fileName to match previous files
qpcr2017_MELT_Lorenzo$fileName <- gsub("_MeltCurve", ".XLS", qpcr2017_MELT_Lorenzo$fileName)

# Manual error : on plate 13.06.2018
qpcr2017_MELT_Lorenzo$Name <- as.character(qpcr2017_MELT_Lorenzo$Name)
qpcr2017_MELT_Lorenzo[qpcr2017_MELT_Lorenzo$fileName %in% "QPCR13.06.2018correct.XLS.csv" &
              qpcr2017_MELT_Lorenzo$Pos %in% paste0("F", 7:12), "Name"] <- "CEWE_AA_0367"

qpcr2017_Lorenzo <- merge(qpcr2017_Lorenzo, 
                 qpcr2017_MELT_Lorenzo[c("fileName", "Name", "Pos", "No..Tm.SYBR", "Tm.x..SYBR")],
                 by = c("fileName", "Name", "Pos"), all = T)

# Remove controls
qpcr2017_Lorenzo <- qpcr2017_Lorenzo[!qpcr2017_Lorenzo$Name %in% c("water", "NTC", "Noiseband","OFF", "Quantification", "SYBR","n/a" ) &
                     !qpcr2017_Lorenzo$Pos %in% c("Threshold level", "Baseline setting"),]

# Remove this plate as all samples were redone, better, later + no melt curve
qpcr2017_Lorenzo <- qpcr2017_Lorenzo[!qpcr2017_Lorenzo$fileName %in% "QPCR04.04.2018.XLS.csv",]

# No melt curves for them:
unique(qpcr2017_Lorenzo[is.na(qpcr2017_Lorenzo$No..Tm.SYBR),"fileName"])

# Remove stupid empty rows
qpcr2017_Lorenzo <- qpcr2017_Lorenzo[!is.na(qpcr2017_Lorenzo$fileName),]

# manual corrections on sample names
qpcr2017_Lorenzo$Name[qpcr2017_Lorenzo$Name %in% "359"] <- "ILWE_AA_0359"
qpcr2017_Lorenzo$Name[qpcr2017_Lorenzo$Name %in% "ILWE_AA_242"] <- "ILWE_AA_0242"
qpcr2017_Lorenzo$Name[qpcr2017_Lorenzo$Name %in% "CEWE_AA_0330#"] <- "CEWE_AA_0330"
qpcr2017_Lorenzo$Name[qpcr2017_Lorenzo$Name %in% "ILWE_AA_0,48"] <- "ILWE_AA_0348"
qpcr2017_Lorenzo$Name[qpcr2017_Lorenzo$Name %in% "ILWE_AA_0134"] <- "ILWE_AA_0379"
qpcr2017_Lorenzo[qpcr2017_Lorenzo$fileName %in% "QPCR07.06.2018part2.XLS.csv" &
          qpcr2017_Lorenzo$Pos %in% paste0("C", 4:12),"Name"] <- gsub("0410", "0409",
                                                             qpcr2017_Lorenzo[qpcr2017_Lorenzo$fileName %in% "QPCR07.06.2018part2.XLS.csv" &
                                                                       qpcr2017_Lorenzo$Pos %in% paste0("C", 4:12),"Name"])
qpcr2017_Lorenzo[qpcr2017_Lorenzo$Name %in% "CEWE_AA_0330#", "Name"] <- "CEWE_AA_0330"
qpcr2017_Lorenzo[qpcr2017_Lorenzo$Name %in% "CEWE_AA_0436*", "Name"] <- "CEWE_AA_0436"
 
qpcr2017_Lorenzo[qpcr2017_Lorenzo$Name %in% "CEWE_AA_0410" & qpcr2017_Lorenzo$Pos %in% c(paste0("D", 7:12)),"Name"] <- "ILWE_AA_0410"

# likely manual mistake
qpcr2017_Lorenzo[qpcr2017_Lorenzo$Pos %in% c("B4", "B5", "B6") &
          qpcr2017_Lorenzo$fileName %in% "QPCR14.06.2018.XLS.csv", "Name"] <- "CEWE_AA_0424"

# Add full name of sample (tissue + mouseID + eimeriaOrmouse primers + plate)
qpcr2017_Lorenzo$fullName <- paste0(qpcr2017_Lorenzo$Name, "_", qpcr2017_Lorenzo$Target.SYBR, "_",  qpcr2017_Lorenzo$fileName)

# Add tissue and Mouse_ID
x <- strsplit(as.character(qpcr2017_Lorenzo$Name), "_", 1)
qpcr2017_Lorenzo$tissue <- sapply( x, "[", 1)
qpcr2017_Lorenzo$Mouse_ID <- paste0("AA_", sapply( x, "[", 3))
rm(x)
##### END Prepare data #####

##### End cleaning ##### 
allSamples <- unique(data.frame(Mouse_ID = qpcr2017_Lorenzo$Mouse_ID,
                         tissue = qpcr2017_Lorenzo$tissue))

### 1. Remove samples with no delta Ct value for mice (means that only one sample worked)
qpcr2017_Lorenzo <- qpcr2017_Lorenzo[-which(is.na(qpcr2017_Lorenzo$Ct.Mean.SYBR) &
                                              qpcr2017_Lorenzo$Target.SYBR == "mouse"),]

# 2. all samples with no CtMean from Eimeria (but with CtMean from mouse, cleaned before) are NEGATIVE
qpcr2017_Lorenzo <- qpcr2017_Lorenzo %>%
  group_by(fullName) %>%
  mutate(status = if_else(is.na(Ct.Mean.SYBR), "negative", "pending"))%>% 
  data.frame()

# 3. Calculate one delta ct per triplicate
qpcr2017_Lorenzo <- qpcr2017_Lorenzo[!duplicated(qpcr2017_Lorenzo[c("fullName")]),]

qpcr2017_Lorenzo <- calculateDeltaCt(qpcr2017_Lorenzo, mergeBy = c("Mouse_ID", "tissue", "Name", "fileName"))

## Average technical replicates
qpcr2017_Lorenzo <- qpcr2017_Lorenzo %>% 
  group_by(Name) %>%   
  summarise(deltaCt_MminusE = mean(deltaCt_MminusE, na.rm = T)) %>% 
  data.frame()

# If NA, means that there was no Eimeria Ct -> negative. Set at -6 baseline
qpcr2017_Lorenzo$deltaCt_MminusE[is.na(qpcr2017_Lorenzo$deltaCt_MminusE)] <- -6

# Final DF
x <- strsplit(as.character(qpcr2017_Lorenzo$Name), "_", 1)
qpcr2017_Lorenzo <- data.frame(Name = qpcr2017_Lorenzo$Name,
                        deltaCt_MminusE = qpcr2017_Lorenzo$deltaCt_MminusE,
                        tissue = sapply( x, "[", 1),
                        Mouse_ID = paste0("AA_", sapply( x, "[", 3)))
rm(x)

qpcr2017_final_Lorenzo <- qpcr2017_Lorenzo["Mouse_ID"]

# Add CEWE
qpcr2017_final_Lorenzo <- merge(qpcr2017_final_Lorenzo,
                                qpcr2017_Lorenzo[
                                  qpcr2017_Lorenzo$tissue %in% "CEWE", 
                                  c("Mouse_ID", "deltaCt_MminusE")],
                                all.x = T)
names(qpcr2017_final_Lorenzo)[names(qpcr2017_final_Lorenzo) %in% "deltaCt_MminusE"] <- "delta_ct_cewe_MminusE"

# Add ILWE
qpcr2017_final_Lorenzo <- merge(qpcr2017_final_Lorenzo,
                                qpcr2017_Lorenzo[
                                  qpcr2017_Lorenzo$tissue %in% "ILWE", 
                                  c("Mouse_ID", "deltaCt_MminusE")],
                           all.x = T)
names(qpcr2017_final_Lorenzo)[names(qpcr2017_final_Lorenzo) %in% "deltaCt_MminusE"] <- "delta_ct_ilwe_MminusE"

qpcr2017_final_Lorenzo <- unique(qpcr2017_final_Lorenzo)

# Add observer
qpcr2017_final_Lorenzo$observer_qpcr <- "Lorenzo" 

# Add year
qpcr2017_final_Lorenzo$year <- 2017

# formate correctly
qpcr2017_final_Lorenzo <- myqPCRformat(qpcr2017_final_Lorenzo)

############ bind full dataset
qpcrData_2016_2017 <- rbind(
  qpcr2016_final_Enas, 
  qpcr2016_final_Mert, 
  qpcr2016_final_Victor, 
  qpcr2016_final_Tabea,
  qpcr2017_final_Lorenzo)

# Correct spaces in names
qpcrData_2016_2017$Mouse_ID <- gsub(" ", "", qpcrData_2016_2017$Mouse_ID)

# Remove duplicates by taking the more recent/reproductible results
doubledSamples <- qpcrData_2016_2017[qpcrData_2016_2017$Mouse_ID %in% 
                                       qpcrData_2016_2017[
                                         duplicated(qpcrData_2016_2017$Mouse_ID), 
                                         "Mouse_ID"], ]
# AA_0064 Tabea
# AA_0066 Tabea
# AA_0070 Victor&Tabea
# AA_0088 Victor&Mert
# AA_0100 Mert
# AA_0139 Mert
toKeep <- rbind(
  doubledSamples[doubledSamples$Mouse_ID %in% "AA_0064" & doubledSamples$observer %in% "Tabea",],
  doubledSamples[doubledSamples$Mouse_ID %in% "AA_0066" & doubledSamples$observer %in% "Tabea",],
  doubledSamples[doubledSamples$Mouse_ID %in% "AA_0100" & doubledSamples$observer %in% "Mert",],
  doubledSamples[doubledSamples$Mouse_ID %in% "AA_0139" & doubledSamples$observer %in% "Mert",])

# ileum stronger
AA_0070 <- cbind(
  doubledSamples[doubledSamples$Mouse_ID %in% "AA_0070" & doubledSamples$observer %in% "Victor", 
                 c("Mouse_ID", "delta_ct_MminusE", "delta_ct_ilwe_MminusE", "qPCRstatus")], 
  delta_ct_cewe_MminusE = doubledSamples[doubledSamples$Mouse_ID %in% "AA_0070" & doubledSamples$observer %in% "Tabea", 
                 c("delta_ct_cewe_MminusE")])
AA_0070$qPCRsummary <- "ileum stronger"
AA_0070$year <- 2016
AA_0070$observer <- "Tabea&Victor"

# ileum stronger
AA_0088 <- cbind(
  doubledSamples[doubledSamples$Mouse_ID %in% "AA_0088" & doubledSamples$observer %in% "Victor", 
                 c("Mouse_ID", "delta_ct_ilwe_MminusE", "qPCRstatus")], 
  doubledSamples[doubledSamples$Mouse_ID %in% "AA_0088" & doubledSamples$observer %in% "Mert", 
                                         c("delta_ct_MminusE", "delta_ct_cewe_MminusE")])
AA_0088$qPCRsummary <- "cecum stronger"
AA_0088$year <- 2016
AA_0088$observer <- "Mert&Victor"

goodOnes <- rbind(toKeep, AA_0070, AA_0088)

# And finally correct
qpcrData_2016_2017 <- qpcrData_2016_2017[-which(qpcrData_2016_2017$Mouse_ID %in% doubledSamples$Mouse_ID),]
qpcrData_2016_2017 <- rbind(qpcrData_2016_2017, goodOnes)

# a little plot to be happy
library(ggplot2)
ggplot(qpcrData_2016_2017, aes(x = Mouse_ID, y = delta_ct_MminusE, 
                               col = qPCRsummary)) +
  geom_point() +
  facet_grid(.~year) +
  theme_bw()

# Write out
write.csv(qpcrData_2016_2017, "../FINALqPCR_2016_2017.csv", row.names = F)
