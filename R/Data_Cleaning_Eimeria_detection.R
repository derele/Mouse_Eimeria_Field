library(ggplot2)

# General data
miceTable <- read.csv("../raw_data/MiceTable_2014to2017.csv")

## Load data from oocysts counting 
data2015 <- read.csv("../raw_data/Eimeria_detection/FINAL2015Oocysts.csv")
data2016 <- read.csv("../raw_data/Eimeria_detection/FINAL2016Oocysts.csv")
data2017 <- read.csv("../raw_data/Eimeria_detection/FINAL2017Oocysts.csv")

myData <- rbind(data2015, data2016, data2017)

# Remove mice without OPG
myData$OPG <- as.numeric(as.character(myData$OPG))
myData <- myData[!is.na(myData$OPG),]

# add HI, GPS coordinates, BCI (and all)
myData <- merge(miceTable, myData, by = "Mouse_ID")

# plot
ggplot(myData) +
  geom_point(aes(x = HI, y = OPG, fill = as.factor(year)),
             pch = 21, alpha = .5, size = 3) +
  theme_bw()

ggplot(myData) +
  geom_point(aes(x = HI, y = log10(OPG)+1, fill = as.factor(year)),
             pch = 21, alpha = .5, size = 3) +
  theme_bw()

# Prevalence / year
prevTable <- as.data.frame.matrix(table(myData$OPG > 0, myData$year))
prevTable[3,] <- round(prevTable[2,] / colSums(prevTable) *100,2)
rownames(prevTable)[3] <- "prevalence(%)"

prevTable

# which ones are juveniles? Can we kick them out for the analysis on Eimeria?
ggplot(myData, aes(x = HI, y = Body_length, col = OPG)) +
  geom_point() +
  scale_color_continuous(low = "grey", high = "red") +
  theme_bw()

myData[myData$Body_length < 60,]

# Density curve : do we catch mice all along the HMHZ?
ggplot(myData, aes(x = HI))+
  geom_density(color = "darkblue", fill = "lightblue") +
  theme_bw()
# Annoying. The model should get this issue
# And what about just the infected ones?
ggplot(myData[myData$OPG>0,], aes(x = HI))+
  geom_density(color = "darkblue", fill = "lightblue") +
  theme_bw()
#same type