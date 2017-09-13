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
merge(counts, counts2, by = "Address", all = TRUE)

## Replace by hand:
correctloc <- function(addresstochange, withnewaddress, codetochange, withnewcode){
  dissection2017$Address[which(dissection2017$Address == addresstochange)] <- withnewaddress
  dissection2017$Address[which(dissection2017$Code == codetochange)] <- withnewcode
}

correctloc(addresstochange, withnewaddress, codetochange, withnewcode)

## Write out new dissection table (beware here if done late :P GIT AT ALL STEPS)
write.csv(x = dissection2017, file = "../raw_data/HZ17_September_Mice_Dissection.csv", row.names = FALSE)

