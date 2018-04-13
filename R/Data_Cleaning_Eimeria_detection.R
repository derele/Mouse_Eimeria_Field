library(ggplot2)

# General data
miceTable <- read.csv("../raw_data/MiceTable_2014to2017.csv")

## Load data from oocysts counting 
data2015 <- read.csv("../raw_data/Eimeria_detection/FINAL2015Oocysts.csv")
data2016 <- read.csv("../raw_data/Eimeria_detection/FINAL2016Oocysts.csv")
data2017 <- read.csv("../raw_data/Eimeria_detection/FINAL2017Oocysts.csv")

myData <- rbind(data2015, data2016, data2017)

# add HI, GPS coordinates, BCI (and all)
myData <- merge(miceTable, myData, by = "Mouse_ID")

# plot
ggplot(myData) +
  geom_point(aes(x = HI, y = log10(as.numeric(as.character(OPG)) +1), fill = as.factor(year)),
             pch = 21, alpha = .5, size = 3) +
  theme_bw() +
  geom_point(data = myData[is.na(as.numeric(as.character(myData$OPG))),],
             aes(x = HI, y = 2), pch = "NA") 

# Prevalence / year
prevTable <- as.data.frame.matrix(table(myData$OPG > 0, myData$year))
prevTable[3,] <- round(prevTable[2,] / colSums(prevTable) *100,2)
rownames(prevTable)[3] <- "prevalence(%)"

prevTable
