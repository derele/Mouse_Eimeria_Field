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

##################### Eimeria detection oocysts flotation ####################
source("../R/functions/addFlotationResults.R")
myData <- addFlotationResults(miceTable)$newDF

# KEEP ONLY ONES WE FLOTATED!!
myData <- myData[!is.na(myData$OPG),]
# MICE NOT FOUND IN miceTable: "SK_3174"
missingHIMice <- c("SK_3174", as.character(myData$Mouse_ID)[is.na(myData$HI)])

prevalenceFlotation <- getPrevalenceTable(myTable = table(myData$OPG > 0,
                                                          myData$year))
prevalenceFlotation

# Create map of samples
mapHMHZ <- HI.map(df = myData, size = 2, alpha = .3, margin = 0.2, zoom = 8) 
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

getPrevalenceTable(table(myData$Ap5_PCR, myData$year))

getPrevalenceTable(table(myData$PCR.positive, myData$year))

#################### Eimeria detection qPCR ####################
source("../R/functions/addqPCRresults.R")
myData <- addqPCRresults(myData)

# the full values are in myData$delta_ct_MminusE
getPrevalenceTable(table(myData$qPCRstatus, myData$year))

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
# 
# ##Plot
# grid.newpage()
# draw.triple.venn(area1 = nrow(subset(finalData, n18S_Seq == "positive")), 
#                  area2 = nrow(subset(finalData, COI_Seq == "positive")), 
#                  area3 = nrow(subset(finalData, ORF470_Seq == "positive")), 
#                  n12 = nrow(subset(finalData, n18S_Seq == "positive"&COI_Seq== "positive")), 
#                  n23 = nrow(subset(finalData, COI_Seq== "positive"&ORF470_Seq =="positive")), 
#                  n13 = nrow(subset(finalData, n18S_Seq == "positive"&ORF470_Seq == "positive")), 
#                  n123 = nrow(subset(finalData, n18S_Seq == "positive"&ORF470_Seq == "positive"&COI_Seq=="positive")),
#                  category = c("18S", "COI", "ORF470"), 
#                  lty = rep(1,3), col = c("dodgerblue4", "firebrick3", "darkgreen"), lwd = rep(2,3),
#                  fill = c("dodgerblue4", "firebrick3", "darkgreen"), alpha = c(0.3, 0.3, 0.3), cex =2, cat.cex = 2.5, cat.default.pos = 'outer', 
#                  cat.col = c("dodgerblue4", "firebrick3", "darkgreen"))
# 
# #Plot
# grid.newpage()
# draw.pairwise.venn(area1= length(which(finalData$Ap5 %in% "positive")), area2= length(which(finalData$Flot %in% "positive")), cross.area = length(which(finalData$Flot %in% "positive" & finalData$Ap5 %in% "positive")))
# 
# 
# #grid.newpage()
# #draw.pairwise.venn(22, 20, 11, category = c("Dog People", "Cat People"), lty = rep("blank", 
# #                                                                                  2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 
# #                                                                                                                                                        0), cat.dist = rep(0.025, 2), scaled = FALSE)
# #https://rstudio-pubs-static.s3.amazonaws.com/13301_6641d73cfac741a59c0a851feb99e98b.html

#################### BCI by HI and OPG #################### 
plotBCI <- ggplot(myData, aes(x = HI, y = BCI, fill = myData$OPG + 1)) +
  geom_point(pch = 21, size = 4, col = "black") +
  scale_fill_gradient(low = "white", high = "black", 
                      breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000), 
                      name = "OPG", trans = "log") +
  theme_bw()

summary(lm(myData$BCI ~ myData$OPG + myData$HI))

