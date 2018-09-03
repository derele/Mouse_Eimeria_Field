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
## Load data from oocysts counting 
flotDF <- read.csv("../raw_data/Eimeria_detection/FINALOocysts2015to2017.csv")
LorenzoDF <- read.csv("../raw_data/Eimeria_detection/Eimeria_oocysts_2015&2017_Lorenzo.csv")
## Import PCR data
PCRdf <- read.csv("../raw_data/Eimeria_detection/Inventory_contents_all.csv")
## Import qPCR data
qpcrData2016 <- read.csv("../raw_data/Eimeria_detection/qPCR_2016.csv")
names(qpcrData2016)[1] <- "Mouse_ID"
qpcrData2017 <- read.csv("../raw_data/Eimeria_detection/qPCR_2017.csv")

#################### Load data ####################
# General data
miceTable <- makeMiceTable("../../Data_important/")

# Remove other rodents
otherRodentsID <- c(miceTable$Mouse_ID[miceTable$Species %in% "Pet mus musculus"],
                    miceTable$Mouse_ID[grep("ZZ", miceTable$Mouse_ID)])
myData <- miceTable[!miceTable$Mouse_ID %in% otherRodentsID,]

##################### Eimeria detection oocysts flotation ####################
flotDF$OPG <- as.numeric(as.character(flotDF$OPG))
flotDF <- flotDF[!is.na(flotDF$OPG),]

## Lorenzo count (in 1mL dilution) for comparison
LorenzoDF <- LorenzoDF[!is.na(LorenzoDF$OPG),]

### Plot comparative Alice (dilution 0.1mL for most samples) and Lorenzo (dilution 1mL)
compData <- merge(flotDF, LorenzoDF, by = "Mouse_ID", all = T)

# merge with current data
myData <- merge(myData, flotDF, all = T)

### Comparison 2 methods of flotation
# How many samples new were detected by decreasing the dilution?
N1 <- sum(compData$OPG.x > 0 & compData$OPG.y == 0, na.rm = T)

adjrsq <- summary(lm(compData$OPG.x ~ compData$OPG.y))$adj.r.squared

plot1 <- ggplot(
  compData, aes(x = OPG.x+1, y = OPG.y+1)) +
  geom_point(alpha = .5, size = 4) +
  coord_equal(ratio=1) +
  xlab("OPG + 1 counted in 0.1ml") + 
  ylab("OPG + 1 counted in 1ml") + 
  geom_smooth(method = "lm", se = FALSE, col = "red") +
  geom_abline(intercept = 0, slope = 1, linetype = 3) +
  scale_y_log10() + 
  scale_x_log10() +
  theme_bw()
plot1

# Remove other rodents
myData <- myData[!myData$Mouse_ID %in% otherRodentsID,]

# correct year
myData$Year[is.na(myData$Year)] <- myData$year[is.na(myData$Year)]
myData <- subset(myData, select = -c(year))

# plot OPG that we keep
plotSmoothOPG <- ggplot(myData[myData$OPG >0,], aes(x = HI, y = OPG+1)) +
  geom_point(aes(fill = as.factor(Year)), pch = 21, alpha = .8, size = 4) +
  geom_smooth(se=F) +
  scale_y_log10() +
  theme_bw() +
  theme(legend.position="top") +
  theme(legend.title = element_blank())
plotSmoothOPG

##################### Eimeria detection PCR ####################

#correct wrong names
toremove <- paste0(paste0("X",1:30, "_"), collapse = "|")
names(PCRdf) <- gsub(toremove, "", names(PCRdf))
names(PCRdf)[names(PCRdf)%in%"ID_mouse"] <- "Mouse_ID"
names(PCRdf)[names(PCRdf)%in%"Year"] <- "yearpcr"
  
PCRdf$Mouse_ID = gsub(" ", "", PCRdf$Mouse_ID) # fix the extra space

# by default, I enter PCRstatus as negative, then overwrite
PCRdf$PCRstatus = "negative"
  
# PCR positive = one of the 3 other markers than AP5 sequenced 
# (Ap5 was used for detection only, the other markers for confirmation)
PCRdf$PCRstatus[PCRdf$`18S_Seq` == "positive" | 
                  PCRdf$COI_Seq == "positive" | 
                  PCRdf$ORF470_Seq == "positive"] <- "positive"

# PCRstatus is NA if everything is NA
PCRdf$PCRstatus[is.na(PCRdf$Ap5_PCR) & 
                  is.na(PCRdf$`18S_Seq`) &
                  is.na(PCRdf$COI_Seq) &
                  is.na(PCRdf$ORF470_Seq)] <- NA

# merge with actual df
myData <- merge(myData, PCRdf, all = T)

# Remove other rodents
myData <- myData[!myData$Mouse_ID %in% otherRodentsID,]

# correct year
myData$Year[is.na(myData$Year)] <- myData$yearpcr[is.na(myData$Year)]
myData <- subset(myData, select = -c(yearpcr))

#################### Eimeria detection qPCR ####################

qpcrData <- qpcrData2016[qpcrData2016$observer_qpcr == "Mert",]

  # # Merge both years
  # # qpcrData <- rbind(qpcrData2016, qpcrData2017Clean)
  # 
  # #####
  # 
  # # Did Enas calculate the other way around? 
  # qpcrData$delta_ct_cewe[qpcrData$observer_qpcr == "Enas"] <- 
  #   - qpcrData$delta_ct_cewe[qpcrData$observer_qpcr == "Enas"]
  # 
  # qpcrData$delta_ct_ilwe[qpcrData$observer_qpcr == "Enas"] <- 
  #   - qpcrData$delta_ct_ilwe[qpcrData$observer_qpcr == "Enas"]
  
# deltaCT = ct eimeria - ct mouse. If high infection, low deltaCT
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

# merge
myData <- merge(myData, qpcrData, by = "Mouse_ID", all = T)
  
# plot qPCR
plotSmoothqPCR <- ggplot(myData[myData$delta_ct_MminusE > - 6,], aes(x = HI)) +
  geom_point(aes(y = delta_ct_MminusE, fill = qPCRsummary), pch = 21, alpha = .8, size = 4) +
  geom_smooth(aes(y = delta_ct_MminusE))+#, col = as.factor(Year))) +
  theme_bw() +
  theme(legend.position="top") +
  theme(legend.title = element_blank()) 
  # facet_grid(Year ~.)
plotSmoothqPCR

# Remark of Justyna Wolinska: some individuals here HAVE qPCR value, but no oocyst count?? Plot

ggplot(myData, aes(x = HI, y = delta_ct_MminusE, col = OPG > 0)) +
  geom_point(size = 3) + 
  geom_hline(yintercept = -6) +
  theme_bw()


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
