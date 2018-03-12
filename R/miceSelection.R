library(ggplot2)

# mouse2016 <- read.csv("../raw_data/HZ16_Mice_18-07-16_dissections.csv") 
mouse2017 <- read.csv("../raw_data/HZ17_September_Mice_Dissection.csv")
alldata <- read.csv("../raw_data/MiceTable_2014to2017.csv")
mouse2015 <- alldata[alldata$Year == 2015,]

miceCountedOo <- read.csv("../raw_data/Eimeria_detection/Alice_newdilution_oocysts_counts_jan2018.csv")
myMice <- miceCountedOo$Mouse_ID[!is.na(miceCountedOo$OPG)]

# Select only counted
mouse2015 <- mouse2015[mouse2015$Mouse_ID %in% myMice,]
mouse2017 <- mouse2017[mouse2017$Mouse_ID %in% myMice,]

# no young
# mouse2016ID <- mouse2016[-grep("young", mouse2016$Sex),"ID_mouse"]
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
# micePicked2016 <- pickMice(mouse2016ID, alldata)
micePicked2017 <- pickMice(mouse2017ID, alldata)

myData <- rbind(micePicked2015, micePicked2017)

# Add OPG info
myData <- merge(myData, data.frame(Mouse_ID = miceCountedOo$Mouse_ID, 
                                   OPG = miceCountedOo$OPG))

# To compare with all
myDataAllInfo <- merge(alldata, data.frame(Mouse_ID = miceCountedOo$Mouse_ID, 
                                           OPG = miceCountedOo$OPG))
myDataAllInfo$location <- paste(round(myDataAllInfo$Longitude,4), round(myDataAllInfo$Latitude, 4))

write.csv(myData, "../raw_data/Eimeria_detection/Partial_mice_usable_for_model.csv", row.names = F)
write.csv(myDataAllInfo, "../raw_data/Eimeria_detection/ALL_mice_usable_for_model.csv", row.names = F)

## Plot and count mice
plotAndInfo <- function(data){
  print("N farms sampled?")
  print(length(unique(data$location)))
  
  print("N mice sampled?")
  print(length(data$Mouse_ID))
  
  print("How many of each sex?")
  print(table(data$Sex))
  
  ggplot(data, aes(x = HI, y = log10(OPG + 1), 
                   col = as.factor(Sex), fill = as.factor(Sex))) +
    geom_point(size = 3) +
    geom_smooth() +
    theme_bw()
}

plotAndInfo(myData)
plotAndInfo(myDataAllInfo)
