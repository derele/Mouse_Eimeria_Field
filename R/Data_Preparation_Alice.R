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

# MICE NOT FOUND IN miceTable: "SK_3174"
missingHIMice <- c(as.character(myData$Mouse_ID[
  !myData$Mouse_ID %in% intersect(miceTable$Mouse_ID, myData$Mouse_ID)]),
  as.character(myData$Mouse_ID)[is.na(myData$HI)])

prevalenceFlotation <- getPrevalenceTable(myTable = table(myData$OPG > 0,
                                                          myData$year))
prevalenceFlotation

# # KEEP ONLY ONES WE FLOTATED!!
datatomap <- myData[!is.na(myData$OPG) & !is.na(myData$Longitude) & !is.na(myData$Latitude),]

# Create map of samples
mapHMHZ <- HI.map(df = datatomap, size = 2, alpha = .3, margin = 0.2, zoom = 8) 
mapHMHZ

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

tabTot <- getPrevalenceTable(table(myData$EimeriaStatus, myData$year))
tabTot

#################### Eimeria detection qPCR ####################
source("../R/functions/addqPCRresults.R")
myData <- addqPCRresults(myData)

# the full values are in myData$delta_ct_MminusE
tabqpcr <- getPrevalenceTable(table(myData$qPCRstatus, myData$year))

#################### General stats on sampling ####################

# mean and 95% ci of N of mice caught / farm (assuming normal distribution)
mean(by(myData, myData["farm"], nrow))
qnorm(0.975)*sd(by(myData, myData["farm"], nrow))/sqrt(nrow(myData))

# density of hybrids
plotDensHI <- ggplot(myData, aes(x = HI)) +
  geom_histogram(binwidth = 0.05, col = "black", fill = "lightblue") +
  theme_bw()
plotDensHI

# Evolution of prevalence / year / farm
prevFarmYear <- ddply(myData,  .(farm, year), function(x){
  return(c(prev = sum(x$OPG > 0) / length(x$OPG) * 100,
           N = length(x$OPG)))})

## Hybrid index calculation:
minHINloci = min(as.numeric(substr(myData$HI_NLoci, 4,6)), na.rm = T)
maxHINloci = max(as.numeric(substr(myData$HI_NLoci, 4,6)), na.rm = T)
meanHINloci = round(mean(as.numeric(substr(myData$HI_NLoci, 4,6)), na.rm = T))

######### Compare our methods of detection
##Venn diagram 

# first, compare PCR and oocysts
completeData1 <- myData[!is.na(myData$OPG) &
                         !is.na(myData$PCRstatus),]
        
myVennDiagram <- function(data){      
  area1 = nrow(subset(completeData1, PCRstatus == "positive"))
  area2 = nrow(subset(completeData1, OPG > 0))
  ## areas of 2-group overlap
  cross.area = nrow(subset(completeData1, PCRstatus == "positive" & OPG > 0))
  
  grid.newpage()
  draw.pairwise.venn(area1 = area1, 
                     area2 = area2, 
                     cross.area = cross.area, 
                     category = c("PCR", "OPG"),
                     fill = c("blue", "red"), 
                     alpha = .5, lwd = 0)
}

myVennDiagram(completeData1)

# to compare, keep only samples tested for the 3 methods
completeData <- myData[!is.na(myData$OPG) &
                         !is.na(myData$PCRstatus) &
                         !is.na(myData$delta_ct_MminusE),]

#################### BCI by HI and OPG #################### 
plotBCI <- ggplot(myData, aes(x = HI, y = BCI, fill = myData$OPG + 1)) +
  geom_point(pch = 21, size = 4, col = "black") +
  scale_fill_gradient(low = "white", high = "black", 
                      breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000), 
                      name = "OPG", trans = "log") +
  theme_bw()

summary(lm(myData$BCI ~ myData$OPG + myData$HI))

