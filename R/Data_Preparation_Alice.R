# Preparation data oocysts count
# Alice Balard
# June 2018
library(ggplot2)
library(ggmap)
library(data.table)
library(plyr)
source("HMHZ_Functions.R")

#################### Load data ####################
# General data
miceTable <- read.csv("../raw_data/MiceTable_2014to2017.csv")

## Load data from oocysts counting 
fullDF <- read.csv("../raw_data/Eimeria_detection/FINALOocysts2015to2017.csv")
fullDF <- fullDF[!is.na(fullDF$OPG),]

## Lorenzo count (in 1mL dilution) for comparison
LorenzoDF <- read.csv("../raw_data/Eimeria_detection/Eimeria_oocysts_2015&2017_Lorenzo.csv")
LorenzoDF <- LorenzoDF[!is.na(LorenzoDF$OPG),]

### Plot comparative Alice (dilution 0.1mL for most samples) and Lorenzo (dilution 1mL)
compData <- merge(fullDF, LorenzoDF, by = "Mouse_ID", all = T)

# read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/raw_data/Eimeria_detection/FINAL2015Oocysts.csv")

# How many samples new were detected by decreasing the dilution?
N1 <- sum(compData$OPG.x > 0 & compData$OPG.y == 0, na.rm = T)

plot1 <- ggplot(compData, aes(x = OPG.x, y = OPG.y)) +
  geom_point(alpha = .5, size = 4) +
  coord_equal(ratio=1) +
  xlab("counted in 0.1ml") + 
  ylab("counted in 1ml") + 
  geom_smooth(method = "lm", se = FALSE, col = "red") +
  geom_abline(intercept = 0, slope = 1, linetype = 3) +
  scale_y_log10() + 
  scale_x_log10() +
  theme_bw()

prevalenceFlotation <- getPrevalenceTable(myTable = table(fullDF$OPG > 0, fullDF$year))
