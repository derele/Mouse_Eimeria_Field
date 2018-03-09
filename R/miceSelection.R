mouse2016 <- read.csv("../raw_data/HZ16_Mice_18-07-16_dissections.csv") 
mouse2017 <- read.csv("../raw_data/HZ17_September_Mice_Dissection.csv")
alldata <- read.csv("../raw_data/MiceTable_2014to2017.csv")
mouse2015 <- alldata[alldata$Year == 2015,]

# no young
mouse2016ID <- mouse2016[-grep("young", mouse2016$Sex),"ID_mouse"]
mouse2017ID <- mouse2017[-grep("young", mouse2017$Status), "Mouse_ID"]
# no info for that for 2015...
mouse2015ID <- mouse2015$Mouse_ID

pickMice <- function(IDs, alldata) {
  # round loc
  subdata <- alldata[alldata$Mouse_ID %in% IDs,]
  subdata$location <- paste(round(subdata$Longitude,4), round(subdata$Latitude, 4))
  # take 1 male, 1 female per location only
  subF <- subdata[subdata$Sex == "F",]
  subM <- subdata[subdata$Sex == "M",]
  subF <- subset(subF, !duplicated(subF$location))
  subM <- subset(subM, !duplicated(subM$location))
  return(rbind(subF, subM))
}

micePicked2015 <- pickMice(mouse2015ID, alldata)
micePicked2016 <- pickMice(mouse2016ID, alldata)
micePicked2017 <- pickMice(mouse2017ID, alldata)

myData <- rbind(micePicked2015, micePicked2016, micePicked2017)
