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
source("HMHZ_Functions.R")

#################### Load data ####################
# General data
miceTable <- read.csv("../raw_data/MiceTable_2014to2017.csv")

# Correct HI error
miceTable$HI[miceTable$HI > 1 & !is.na(miceTable$HI)] <- 
  miceTable$HI[miceTable$HI > 1 & !is.na(miceTable$HI)]/1000

# Correct Lat/Lon errors
miceTable$Longitude[miceTable$Longitude > 100 & !is.na(miceTable$Longitude)] <- 
  miceTable$Longitude[miceTable$Longitude > 100 & !is.na(miceTable$Longitude)] / 1000

miceTable$Latitude[miceTable$Latitude > 100 & !is.na(miceTable$Latitude)] <- 
  miceTable$Latitude[miceTable$Latitude > 100 & !is.na(miceTable$Latitude)] / 1000

# Cluster by localities: rounded to about 700 meters ???...
# miceTable$Latitude <- round(miceTable$Latitude, 2)
# miceTable$Longitude <- round(miceTable$Longitude, 2)

# Body weight correction of x 1000 error
miceTable$Body_weight[!is.na(miceTable$Body_weight) &
                        miceTable$Body_weight > 100] <-
  miceTable$Body_weight[!is.na(miceTable$Body_weight) &
                          miceTable$Body_weight > 100] / 1000

# Body condition index as log body mass/log body length (Hayes et al. 2014)
miceTable$BCI <- log(miceTable$Body_weight) / log(miceTable$Body_length)

# add farm (TODO better localisation)
miceTable$farm <- paste0(miceTable$Longitude, miceTable$Latitude)

## remove empty rows
miceTable <- miceTable[!is.na(miceTable$Mouse_ID),]

##################### Eimeria detection oocysts flotation ####################
source("../R/functions/addFlotationResults.R")
myData <- addFlotationResults(miceTable)$newDF

# correct year
myData$year[is.na(myData$year)] <- myData$Year[is.na(myData$year)]
myData <- subset(myData, select = -c(Year))
myData$Mouse_ID[is.na(myData$year)] # check, must be null

# prevalenceFlotation <- getPrevalenceTable(myTable = table(myData$OPG > 0,
#                                                           myData$year))
# prevalenceFlotation

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
  geom_point(aes(fill = as.factor(year)), pch = 21, alpha = .8, size = 4) +
  geom_smooth(se=F) +
  scale_y_log10() +
  theme_bw() +
  theme(legend.position="top") +
  theme(legend.title = element_blank())
plotSmoothOPG

##################### Eimeria detection PCR ####################
source("../R/functions/addPCRresults.R")
myData <- addPCRresults(myData)

# correct year
myData$year <- myData$year.x
myData$year[is.na(myData$year)] <- myData$year.y[is.na(myData$year)]
myData <- subset(myData, select = -c(year.x, year.y))
myData$Mouse_ID[is.na(myData$year)] # check, must be null

# JUST with Ap5
# getPrevalenceTable(table(myData$Ap5_PCR, myData$year))

# Flotation positive OR PCR positive (= full Eimeria status)
# Set to NA
myData$EimeriaStatus <- NA
# Negative if one test is
myData$EimeriaStatus[myData$PCRstatus == "negative" | myData$OPG <= 0] <- "negative"
# But overwrite positive if one of the test is :)
myData$EimeriaStatus[myData$PCRstatus == "positive" | myData$OPG > 0] <- "positive"

# tabTot <- getPrevalenceTable(table(myData$EimeriaStatus, myData$year))
# tabTot

#################### Eimeria detection qPCR ####################
source("../R/functions/addqPCRresults.R")
myData <- addqPCRresults(myData)

# # the full values are in myData$delta_ct_MminusE
# tabqpcr <- getPrevalenceTable(table(myData$qPCRstatus, myData$year))
# tabqpcr

########### TO CORRECT ########### 
## Which mice are not found?

# 1. no HI given for these mice
missingHIMice <- myData$Mouse_ID[is.na(myData$HI)]

# wildpark Schorfheide (not needed, test)
wsh <- c(paste0("AA_000", 1:9), paste0("AA_00", 10:46))
# apodemus caught in 2016
apd <- c("A_0001", "A_0002", "A_0003")
# useless info
useless <- c(wsh, apd)

# 2. these mice were not in the initial data from jarda
notfound <- myData$Mouse_ID[!myData$Mouse_ID %in% miceTable$Mouse_ID]

# 3. total mice missing
miceInfoNeeded <- unique(as.character(notfound), as.character(missingHIMice))[
  !unique(as.character(notfound), as.character(missingHIMice)) %in% useless]

#################### General stats on sampling ####################

# Keep mice with OPG, PCR or qPCR status, from the Brandenburg-MVP transect
myDataStudyAlice <- myData[!myData$Mouse_ID %in% useless & # no wilpark schofheide
                             (!is.na(myData$OPG) | 
                                !is.na(myData$PCRstatus) |
                                !is.na(myData$qPCRstatus)),] # have one of these status

# latitude or longitude missing for mice:
latLongMissing <- myDataStudyAlice$Mouse_ID[
  is.na(myDataStudyAlice$Latitude) |
    is.na(myDataStudyAlice$Longitude)]

# keep only North Germany
myDataStudyAlice <- myDataStudyAlice[!is.na(myDataStudyAlice$Latitude) &
                                       myDataStudyAlice$Latitude > 51, ]

# Total
Nmice <- nrow(myDataStudyAlice)
Nfarm <- length(unique(myDataStudyAlice$farm))

# Create map of samples
mapHMHZ <- HI.map(df = myDataStudyAlice, size = 2, alpha = .3, margin = 0.2, zoom = 8) 
mapHMHZ

# mean and 95% ci of N of mice caught / farm (assuming normal distribution)
MEAN <- mean(by(myDataStudyAlice, myDataStudyAlice["farm"], nrow))
CI <- qnorm(0.975)*sd(by(myDataStudyAlice, myDataStudyAlice["farm"], nrow))/
  sqrt(nrow(myDataStudyAlice))

# density of hybrids
plotDensHI <- ggplot(myDataStudyAlice, aes(x = HI)) +
  geom_histogram(binwidth = 0.05, col = "black", fill = "lightblue") +
  theme_bw()
plotDensHI

# Hybrid index calculation:
minHINloci = min(as.numeric(substr(myDataStudyAlice$HI_NLoci, 4,6)), na.rm = T)
maxHINloci = max(as.numeric(substr(myDataStudyAlice$HI_NLoci, 4,6)), na.rm = T)
meanHINloci = round(mean(as.numeric(substr(myDataStudyAlice$HI_NLoci, 4,6)), na.rm = T))

#Prevalence compared
prevalenceFlotation <- getPrevalenceTable(table(myDataStudyAlice$OPG > 0, myDataStudyAlice$year))
prevalencePCR <- getPrevalenceTable(table(myDataStudyAlice$PCRstatus, myDataStudyAlice$year))
prevalenceqPCR <- getPrevalenceTable(table(myDataStudyAlice$qPCRstatus, myDataStudyAlice$year))

######### Compare our methods of detection ######### 
noo <- table(!is.na(myDataStudyAlice$OPG))[2]
npcr <- table(!is.na(myDataStudyAlice$PCRstatus))[2]
nqPCR <- table(!is.na(myDataStudyAlice$qPCRstatus))[2]

# Positive for qPCR but nothing else?
myDataStudyAlice[which(myDataStudyAlice$qPCRstatus == "positive" & 
                            !is.na(myDataStudyAlice$qPCRstatus) &
                            myDataStudyAlice$PCRstatus == "negative" &
                            myDataStudyAlice$OPG == 0), ]
##Venn diagram 

# first, compare PCR and oocysts
completeData1 <- myDataStudyAlice[!is.na(myDataStudyAlice$OPG) &
                                    !is.na(myDataStudyAlice$PCRstatus),]

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
ggplot(myDataStudyAlice[!is.na(myDataStudyAlice$qPCRstatus),],
       aes(x = delta_ct_MminusE, y = OPG)) +
  geom_point(aes(fill = qPCRsummary),
             alpha = .5, size = 4, pch = 21) +
  geom_smooth(method = "lm", se = FALSE, col = "red") +
  theme_bw()

plotcompOPGqPCR <- ggplot(myDataStudyAlice[!is.na(myDataStudyAlice$delta_ct_MminusE) &
                                             myDataStudyAlice$OPG > 0 &
                                             myDataStudyAlice$delta_ct_MminusE > -6,],
                          aes(x = delta_ct_MminusE, y = OPG)) +
  geom_point(aes(fill = qPCRsummary),
             alpha = .5, size = 4, pch = 21) +
  geom_smooth(method = "lm", se = FALSE, col = "red") +
  scale_y_log10() +
  theme_bw()

data1 <- myDataStudyAlice[!is.na(myDataStudyAlice$delta_ct_MminusE) &
                            !is.na(myDataStudyAlice$OPG), ]
summary(lm(data1$OPG ~ data1$delta_ct_MminusE))

data2 <- data1[data1$OPG > 0 & data1$delta_ct_MminusE > -6 , ]
summary(lm(data2$OPG ~ data2$delta_ct_MminusE))

# to compare, keep only samples tested for the 3 methods
completeData <- myDataStudyAlice[!is.na(myDataStudyAlice$OPG) &
                                   !is.na(myDataStudyAlice$PCRstatus) &
                                   !is.na(myDataStudyAlice$qPCRstatus),]

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
