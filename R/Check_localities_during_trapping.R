## Alice September 2017

## Import files to cross check
dissection2017 <- read.csv("../raw_data/HZ17_September_Mice_Dissection.csv")
dissection2017$Address <- as.character(dissection2017$Address)

trapping2017 <- read.csv("../raw_data/HZ17_Mice_Trap.csv")

## Dissection localities frequency table:
library(plyr)
counts <- ddply(dissection2017, .(dissection2017$Address, dissection2017$Code,
                                  dissection2017$Latitude, dissection2017$Longitude), nrow)
names(counts) <- c("Address", "Code", "Latitude", "Longitude")
counts

## Trapping localities frequency table:
counts2 <- ddply(trapping2017, .(trapping2017$Address, trapping2017$Code,
                                  trapping2017$Latitude, trapping2017$Longitude), nrow)
names(counts2) <- c("Address", "Code", "Latitude", "Longitude")
counts2

## Merge for a visual check:
merge(counts, counts2, by = c("Longitude", "Latitude"), all = TRUE)

## Replace by hand:
correctloc <- function(addresstochange, withnewaddress, codetochange, withnewcode){
  dissection2017$Address[which(dissection2017$Address == addresstochange)] <- withnewaddress
  dissection2017$Address[which(dissection2017$Code == codetochange)] <- withnewcode
}

as.character(dissection2017$Code)

correctloc(addresstochange, withnewaddress, codetochange, withnewcode)

###############################
## Part 2: draw a map
dissection2017$OldOrNew <- "OldLoc"
dissection2017[which(dissection2017$Code == "NA"), ]#$OldOrNew <- "NewLoc"

dissection2017$musornot <- "Mus musculus"
dissection2017[which(dissection2017$Species != "Mus musculus"), ]$musornot <- "other species"

library(ggmap)

dissection2017$Latitude <- as.numeric(dissection2017$Latitude)
dissection2017$Longitude <- as.numeric(dissection2017$Longitude)

# get a map
margin <- 0.2

area <- get_map(location =
                  c(min(na.omit(dissection2017$Longitude) - margin),
                    min(na.omit(dissection2017$Latitude) - margin),
                    max(na.omit(dissection2017$Longitude) + margin),
                    max(na.omit(dissection2017$Latitude) + margin)),
                source = "stamen", maptype="toner-lite",
                zoom = 9)

#plot the map :
ggmap(area) +
  geom_jitter(data = dissection2017, shape = 21, size = 3, 
             aes(Longitude, Latitude, col = musornot), width =  0.005, height = 0.005)
