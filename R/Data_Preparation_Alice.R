# Preparation data oocysts count
# Alice Balard
# June 2018
library(ggplot2)
library(ggmap)
library(data.table)
library(plyr)
library("VennDiagram")
library(grid)
library(gridExtra)
source("functions/HMHZ_Functions.R")
source("functions/makeMiceTable.R")

#################### Load data ####################
# General data
miceTable <- makeMiceTable("../../Data_important/")

# Remove other rodents
otherRodentsID <- c(miceTable$Mouse_ID[miceTable$Species %in% "Pet mus musculus"],
                    miceTable$Mouse_ID[grep("ZZ", miceTable$Mouse_ID)])
myData <- miceTable[!miceTable$Mouse_ID %in% otherRodentsID,]

##################### Eimeria detection oocysts flotation ####################
source("../R/functions/addFlotationResults.R")
myData <- addFlotationResults(myData)$newDF

# Remove other rodents
myData <- myData[!myData$Mouse_ID %in% otherRodentsID,]

# correct year
myData$Year[is.na(myData$Year)] <- myData$year[is.na(myData$Year)]
myData <- subset(myData, select = -c(year))

# How many new samples were found given the new dilution (0.1mL)
N1 <- comparisonFlot(myData)$N1
N1

# What is the adjusted R square between these 2 dilutions
adjrsq <- comparisonFlot(myData)$adjrsq
adjrsq

# Plot the comparison
plot1 <- plotCompData(myData)
plot1

# plot OPG that we keep
plotSmoothOPG <- ggplot(myData[myData$OPG >0,], aes(x = HI, y = OPG+1)) +
  geom_point(aes(fill = as.factor(Year)), pch = 21, alpha = .8, size = 4) +
  geom_smooth(se=F) +
  scale_y_log10() +
  theme_bw() +
  theme(legend.position="top") +
  theme(legend.title = element_blank())
plotSmoothOPG

plotSmoothOPGall <- ggplot(myData[myData$Year %in% 2015:2017,], 
                           aes(x = HI, y = OPG+1)) +
  geom_point(aes(fill = Year), pch = 21, alpha = .8, size = 4) +
  geom_smooth(col = "black") +
  scale_y_log10() +
  theme_bw() +
  theme(legend.position="top") +
  theme(legend.title = element_blank())
plotSmoothOPGall

##################### Eimeria detection PCR ####################
source("../R/functions/addPCRresults.R")
myData <- addPCRresults(myData)

# Remove other rodents
myData <- myData[!myData$Mouse_ID %in% otherRodentsID,]

# correct year
myData$Year[is.na(myData$Year)] <- myData$year[is.na(myData$Year)]
myData <- subset(myData, select = -c(year))

#################### Eimeria detection qPCR ####################
source("../R/functions/addqPCRresults.R")
myData <- addqPCRresults(myData, 
                         pathtoqPCR2016 = "../raw_data/Eimeria_detection/qPCR_2016.csv",
                         pathtoqPCR2017 = "../raw_data/Eimeria_detection/qPCR_2017.csv")

# Remove other rodents
myData <- myData[!myData$Mouse_ID %in% otherRodentsID,]

# plot qPCR
plotSmoothqPCR <- ggplot(myData[myData$delta_ct_cewe < 6,], aes(x = HI)) +
  geom_point(aes(y = -delta_ct_cewe, fill = as.factor(Year)), pch = 21, alpha = .8, size = 4) +
  geom_smooth(aes(y = -delta_ct_cewe))+#, col = as.factor(Year))) +
  theme_bw() +
  theme(legend.position="top") +
  theme(legend.title = element_blank())
plotSmoothqPCR

# Fit the model!! the smooth is not adapted

#################### General stats on sampling ####################
## Remove useless mice:

# wildpark Schorfheide (not needed, test)
wsh <- c(paste0("AA_000", 1:9), paste0("AA_00", 10:46))
# apodemus caught in 2016
apd <- c("A_0001", "A_0002", "A_0003")
# useless info
useless <- c(wsh, apd)

# Keep mice with OPG, PCR or qPCR status
myData <- myData[!myData$Mouse_ID %in% useless & # no wilpark schofheide
                             (!is.na(myData$OPG) | 
                                !is.na(myData$PCRstatus) |
                                !is.na(myData$qPCRstatus)),] # have one of these status

## Which mice are not found (no HI given for these mice)?
miceInfoNeeded <- myData$Mouse_ID[is.na(myData$HI)]

# latitude or longitude missing for mice:
latLongMissing <- myData$Mouse_ID[
  is.na(myData$Latitude) |
    is.na(myData$Longitude)]

# keep only North Germany
myData <- myData[!is.na(myData$Latitude) &
                                       myData$Latitude > 51 &
                                       myData$Longitude < 17, ]

# Total
Nmice <- nrow(myData)
Nfarm <- length(unique(myData$farm))

# Create map of samples
mapHMHZ <- HI.map(df = myData, size = 2, alpha = .3, margin = 0.2, zoom = 8) 
mapHMHZ

# mean and 95% ci of N of mice caught / farm (assuming normal distribution)
MEAN <- mean(by(myData, myData["farm"], nrow))
CI <- qnorm(0.975)*sd(by(myData, myData["farm"], nrow))/
  sqrt(nrow(myData))

# density of hybrids
plotDensHI <- ggplot(myData, aes(x = HI)) +
  geom_histogram(binwidth = 0.05, col = "black", fill = "lightblue") +
  theme_bw()
plotDensHI

# Hybrid index calculation:
minHINloci = min(as.numeric(substr(myData$HI_NLoci, 4,6)), na.rm = T)
maxHINloci = max(as.numeric(substr(myData$HI_NLoci, 4,6)), na.rm = T)
meanHINloci = round(mean(as.numeric(substr(myData$HI_NLoci, 4,6)), na.rm = T))

#Prevalence compared
prevalenceFlotation <- getPrevalenceTable(table(myData$OPG > 0, 
                                                myData$Year))
prevalencePCR <- getPrevalenceTable(table(myData$PCRstatus, 
                                          myData$Year))
prevalenceqPCR <- getPrevalenceTable(table(myData$qPCRstatus, 
                                           myData$Year))

myData$allDetectionMethod <- NA
myData$allDetectionMethod[myData$OPG >0 |
                            myData$PCRstatus == "positive" |
                            myData$qPCRstatus == "positive"] <- "positive"
myData$allDetectionMethod[myData$OPG == 0 &
                            myData$PCRstatus == "negative" &
                            myData$qPCRstatus == "negative"] <- "negative"

prevalenceTot <- getPrevalenceTable(table(myData$allDetectionMethod,
                                            myData$Year ))



######### Compare our methods of detection ######### 
noo <- table(!is.na(myData$OPG))[2]
npcr <- table(!is.na(myData$PCRstatus))[2]
nqPCR <- table(!is.na(myData$qPCRstatus))[2]

# Positive for qPCR but nothing else?
myData[which(myData$qPCRstatus == "positive" & 
                            !is.na(myData$qPCRstatus) &
                            myData$PCRstatus == "negative" &
                            myData$OPG == 0), ]
##Venn diagram 

# first, compare PCR and oocysts
completeData1 <- myData[!is.na(myData$OPG) &
                                    !is.na(myData$PCRstatus),]

myVennDiagram2 <- function(data){      
  area1 = nrow(subset(data, PCRstatus == "positive"))
  area2 = nrow(subset(data, OPG > 0))
  ## areas of 2-group overlap
  cross.area = nrow(subset(data, PCRstatus == "positive" & OPG > 0))
  grid.newpage()
  draw.pairwise.venn(area1 = area1, 
                     area2 = area2, 
                     cross.area = cross.area, 
                     category = c("PCR", "OPG"),
                     col = "transparent", 
                     fill = c("grey","green"),
                     alpha = 0.50,                  
                     cex = 1.5, cat.cex = 1.5, fontfamily = "serif", fontface = "bold",
                     cat.col = c("grey","green"),
                     cat.fontfamily = "serif")
}

myVennDiagram2(completeData1)

# Compare qPCR results and OPG
ggplot(myData[!is.na(myData$qPCRstatus),],
       aes(x = delta_ct_MminusE, y = OPG)) +
  geom_point(aes(fill = qPCRsummary),
             alpha = .5, size = 4, pch = 21) +
  geom_smooth(method = "lm", se = FALSE, col = "red") +
  theme_bw()

plotcompOPGqPCR <- ggplot(myData[!is.na(myData$delta_ct_MminusE) &
                                             myData$OPG > 0 &
                                             myData$delta_ct_MminusE > -6,],
                          aes(x = delta_ct_MminusE, y = OPG)) +
  geom_point(aes(fill = qPCRsummary),
             alpha = .5, size = 4, pch = 21) +
  geom_smooth(method = "lm", se = FALSE, col = "red") +
  scale_y_log10() +
  theme_bw()

data1 <- myData[!is.na(myData$delta_ct_MminusE) &
                            !is.na(myData$OPG), ]
summary(lm(data1$OPG ~ data1$delta_ct_MminusE))

data2 <- data1[data1$OPG > 0 & data1$delta_ct_MminusE > -6 , ]
summary(lm(data2$OPG ~ data2$delta_ct_MminusE))

# to compare, keep only samples tested for the 3 methods
completeData <- myData[!is.na(myData$OPG) &
                                   !is.na(myData$PCRstatus) &
                                   !is.na(myData$qPCRstatus),]

myVennDiagram3 <- function(data){      
  area1 = nrow(subset(data, PCRstatus == "positive"))
  area2 = nrow(subset(data, OPG > 0))
  area3 = nrow(subset(data, qPCRstatus == "positive"))
  ## areas of 2-group overlap
  n12 = nrow(subset(data, PCRstatus == "positive" & OPG > 0)) 
  n23 = nrow(subset(data, qPCRstatus == "positive" & OPG > 0))
  n13 = nrow(subset(data, PCRstatus == "positive" & qPCRstatus == "positive"))
  ## areas of 3-group overlap
  n123 = nrow(subset(data, PCRstatus == "positive" & 
                       OPG > 0 &
                       qPCRstatus == "positive"))
  grid.newpage()
  draw.triple.venn(area1 = area1, area2 = area2, area3 = area3,
                   n12 = n12, n23 = n23, n13 = n13,
                   n123 = n123,
                   category = c("PCR", "OPG", "qPCR"),
                   col = "transparent", 
                   fill = c("grey","green","orange"),
                   alpha = 0.50,                  
                   cex = 1.5, cat.cex = 1.5, fontfamily = "serif", fontface = "bold",
                   cat.col = c("grey","green","orange"),
                   cat.fontfamily = "serif")
}

myVennDiagram3(completeData)

#################### BCI by HI and OPG #################### 
plotBCI <- ggplot(myData, aes(x = HI, y = BCI, fill = myData$OPG + 1)) +
  geom_point(pch = 21, size = 4, col = "black") +
  scale_fill_gradient(low = "white", high = "black", 
                      breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000), 
                      name = "OPG", trans = "log") +
  theme_bw()

summary(lm(myData$BCI ~ myData$OPG + myData$HI))

#### Bonus: mice SNP chip
# infectedMice
myDataHI <- myData[!is.na(myData$HI),]

# positive for flotation and have an hybrid index
N_OPGPositive <- nrow(myDataHI[myDataHI$OPG > 0 & 
                                 !is.na(myDataHI$OPG),])
totalOPG <- nrow(myDataHI[!is.na(myDataHI$OPG),])

# positive for flotation and have an hybrid index
N_qPCRPositive <- nrow(myDataHI[myDataHI$qPCRstatus == "positive" &
                                  !is.na(myDataHI$qPCRstatus), ])
totalqPCR <- nrow(myDataHI[!is.na(myDataHI$qPCRstatus),])

# positive for flotation and have an hybrid index
N_QuantitativePositive <- nrow(myDataHI[(myDataHI$OPG > 0 | 
                                           myDataHI$qPCRstatus == "positive") &
                                          (!is.na(myDataHI$qPCRstatus) |
                                          !is.na(myDataHI$OPG)),])
totalQuantitative <- nrow(myDataHI[!is.na(myDataHI$OPG) | 
                                     !is.na(myDataHI$qPCRstatus),])
# 
# * `r N_QuantitativePositive` out of `r totalQuantitative` are positive for either flotation or qPCR, and have an hybrid index.

# to fix here!!
