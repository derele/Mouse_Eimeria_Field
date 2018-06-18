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

# Correct Lat/Lon errors
miceTable$Longitude[miceTable$Longitude > 100 & !is.na(miceTable$Longitude)] <- 
  miceTable$Longitude[miceTable$Longitude > 100 & !is.na(miceTable$Longitude)] / 1000

miceTable$Latitude[miceTable$Latitude > 100 & !is.na(miceTable$Latitude)] <- 
  miceTable$Latitude[miceTable$Latitude > 100 & !is.na(miceTable$Latitude)] / 1000

# Cluster by localities: rounded to about 700 meters ???...
# miceTable$Latitude <- round(miceTable$Latitude, 2)
# miceTable$Longitude <- round(miceTable$Longitude, 2)

# Create map of samples
# HI.map(df = miceTable, size = 4, alpha = 1, margin = 0.2) 
# 
# HI.map <- function(df, size = 3, margin = 2, zoom = 7, alpha = 0.5, 
#                    source = "stamen", maptype = "toner-lite"){
# #   # get a map
# margin = 0.2
# zoom = 7
# df = miceTable
# https://github.com/dkahle/ggmap/issues/107
# area <- get_map(
#                 source = "stamen", maptype = "toner-lite", zoom = 10)
# 
# 
# get_stamenmap(bbox = 
#                 c(min(df$Longitude[!is.na(df$Longitude)] - margin),
#                   min(df$Latitude[!is.na(df$Latitude)] - margin),
#                   max(df$Longitude[!is.na(df$Longitude)] + margin),
#                   max(df$Latitude[!is.na(df$Latitude)] + margin)), 
#               maptype = "terrain")
# 
# #plot the map :
#   ggmap(area) +
#     geom_point(data = df, shape = 21, size = size,
#                aes(Longitude, Latitude, fill = HI), alpha = alpha) + # set up the points
#     scale_fill_gradient("Hybrid\nindex", high="red",low="blue")   # set up the HI colors
# }  

# Body weight correction of x 1000 error
miceTable$Body_weight[!is.na(miceTable$Body_weight) &
                        miceTable$Body_weight > 100] <-
  miceTable$Body_weight[!is.na(miceTable$Body_weight) &
                          miceTable$Body_weight > 100] / 1000

# Body condition index as log body mass/log body length (Hayes et al. 2014)
miceTable$BCI <- log(miceTable$Body_weight) / log(miceTable$Body_length)

## Load data from oocysts counting 
fullDF <- read.csv("../raw_data/Eimeria_detection/FINALOocysts2015to2017.csv")
fullDF <- fullDF[!is.na(fullDF$OPG),]

## Lorenzo count (in 1mL dilution) for comparison
LorenzoDF <- read.csv("../raw_data/Eimeria_detection/Eimeria_oocysts_2015&2017_Lorenzo.csv")
LorenzoDF <- LorenzoDF[!is.na(LorenzoDF$OPG),]

### Plot comparative Alice (dilution 0.1mL for most samples) and Lorenzo (dilution 1mL)
compData <- merge(fullDF, LorenzoDF, by = "Mouse_ID", all = T)

# How many samples new were detected by decreasing the dilution?
N1 <- sum(compData$OPG.x > 0 & compData$OPG.y == 0, na.rm = T)
N1

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
plot1

adjrsq <- summary(lm(formula = compData$OPG.x ~ compData$OPG.y))$adj.r.squared

prevalenceFlotation <- getPrevalenceTable(myTable = table(fullDF$OPG > 0, fullDF$year))
prevalenceFlotation

# Add HI, GPS coordinates, BCI (and all)
myData <- merge(miceTable, fullDF, by = "Mouse_ID")

# Add Eimeria infection status TO BE DISCUSSED !!!
# myData <- merge(Eimeria[c("Mouse_ID", "PCRpos3markers")], MiceTable, all = T)

# MICE NOT FOUND IN miceTable: "SK_3174"
missingHIMice <- c("SK_3174", as.character(myData$Mouse_ID)[is.na(myData$HI)])

# add farm (TODO better localisation)
myData$farm <- paste0(myData$Longitude, myData$Latitude)

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

#################### Eimeria detection OPG ####################
# plot
plotSmoothOPG <- ggplot(myData[myData$OPG >0,], aes(x = HI, y = OPG+1)) +
  geom_point(aes(fill = as.factor(year)), pch = 21, alpha = .8, size = 4) +
  geom_smooth(se=F) +
  scale_y_log10() +
  theme_bw() +
  theme(legend.position="top") +
  theme(legend.title = element_blank())
plotSmoothOPG

##################### Eimeria detection PCR ####################

# Eimeria <- read.csv("../raw_data/Eimeria_detection/Summary_eimeria.csv")

#################### Eimeria detection qPCR ####################

#################### BCI by HI and OPG #################### 
plotBCI <- ggplot(myData, aes(x = HI, y = BCI, fill = myData$OPG + 1)) +
  geom_point(pch = 21, size = 4, col = "black") +
  scale_fill_gradient(low = "white", high = "black", 
                      breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000), 
                      name = "OPG", trans = "log") +
  theme_bw()

summary(lm(myData$BCI ~ myData$OPG + myData$HI))

