lolo <- read.csv("../LorenzoRAW/CSVFiles/TotalLorenzo.csv", stringsAsFactors = F,
                 na.strings = c("NA", "", " ", "-"))

# Column file name
names(lolo)[names(lolo) == "QPCR01.06.2018.XLS.csv"] <- "fileName"
lolo[,names(lolo) %in% c("fileName", "Target.SYBR")] <- 
  gsub(" ", "", lolo[,names(lolo) %in% c("fileName", "Target.SYBR")] )

# Remove samples with no Ct value
lolo <- lolo[!is.na(lolo$Ct.SYBR),]

# Remove inside headers
lolo <- lolo[lolo$Name != "Name",]

# Remove controls (were used before)
lolo <- lolo[!lolo$Name %in% c("water", "NTC"),]

# manual correction
lolo$Name[lolo$Name == "359"] <- "ILWE_AA_0359"
lolo$Name[lolo$Name == "ILWE_AA_242"] <- "ILWE_AA_0242"
lolo$Name[lolo$Name == "ILWE_AA_0,48"] <- "ILWE_AA_0348"
lolo$Name[lolo$Name == "ILWE_AA_0134"] <- "ILWE_AA_0379"

lolo[lolo$Pos %in% c("B4", "B5", "B6") & 
       lolo$fileName == "QPCR14.06.2018.XLS.csv", "Name"] <- "CEWE_AA_0424"

levels(lolo$fileName)
# Add tissue and EH_ID
x <-strsplit(as.character(lolo$Name), "_", 1)

lolo$tissue <- sapply( x, "[", 1)

lolo$EH_ID <- paste0("AA_", sapply( x, "[", 3))

table(lolo$fileName) #QPCR03.05.2018.XLS.csv not complete!!! cf key usb Victor gave me

# Count triplicates
library(plyr)
summary <- plyr::count(lolo[, c("EH_ID", "tissue", "Target.SYBR", "fileName")])

toCheck <- plyr::count(lolo[, c("EH_ID", "tissue", "Target.SYBR")])

# find the samples with too many repeats and look manually
tooManyRepeats <- toCheck[toCheck$freq > 3,]

alorsOnmerge <- merge(tooManyRepeats, lolo, all.x = T)

alorsOnmerge <- alorsOnmerge[c("EH_ID", "tissue", "Target.SYBR", "fileName")]

alorsOnmerge2 <- unique(alorsOnmerge)

write.csv(alorsOnmerge2, "files to check", row.names = F)

# in which files?
unique(alorsOnmerge2$fileName)

################ After Correction #####################
okRepeats <- toCheck[toCheck$freq <= 3,]
incompleteData <- merge(lolo, okRepeats)

# Remove triplicates (MeanCt was calculated by the program)
myDF <- incompleteData[c("Target.SYBR", "tissue", "EH_ID", "Ct.Mean.SYBR")]
myDF <- unique(myDF)

# remove attempts with no mean (< 2 success. See if we take that or < 3myDF
myDF <- myDF[!is.na(myDF$Ct.Mean.SYBR),]

# check
table(data.frame(table(myDF$EH_ID, myDF$tissue, myDF$Target.SYBR))["Freq"] )

# Calculate deltaCt (Eimeria - Mouse)
for (i in 1:nrow(myDF)){
  if (myDF$Target.SYBR[i] == "mouse"){
    myDF$Ct.Mean.SYBR.Mus[i] <- myDF$Ct.Mean.SYBR[i]
  } else {
    myDF$Ct.Mean.SYBR.Eimeria[i] <- myDF$Ct.Mean.SYBR[i]
  }
}


as.character(myDF$Target.SYBR)

##################### clean

table(lolo$EH_ID, lolo$tissue, lolo$Target.SYBR)

length(unique(lolo$EH_ID, lolo$tissue, lolo$Target.SYBR))
length(unique(lolo$EH_ID))


table(lolo$EH_ID, lolo$fileName)



lolo$Ct.Mean.SYBR <- as.numeric(as.character(lolo$Ct.Mean.SYBR))


library(ggplot2)
ggplot(lolo, aes(x = lolo$EH_ID, y = lolo$Ct.Mean.SYBR)) +
  geom_point(aes(col = lolo$Target.SYBR))

table(lolo$EH_ID, lolo$tissue, lolo$Target.SYBR)

write.csv(lolo, file = "/home/alice/Schreibtisch/git_projects/Mouse_Eimeria_Databasing/raw_data/Eimeria_detection/raw_qPCR/lorenzo2017.csv", row.names = F)

library(dplyr)

newLolo <- lolo %>% 
  group_by(Name, Target.SYBR, Ct.Mean.SYBR) %>% 
  mutate(Ct.SYBR_Rep1 = Ct.SYBR[1],
         Ct.SYBR_Rep2 = Ct.SYBR[2],
         Ct.SYBR_Rep3 = Ct.SYBR[3]) %>%
  data.frame()

newLolo <- newLolo[
  as.numeric(
    row.names(
      unique(newLolo[c("Name", "Target.SYBR", "Ct.Mean.SYBR")]))),]

newLolo <- newLolo[!names(newLolo) %in% "Ct.SYBR"]

table(newLolo$EH_ID, newLolo$tissue, newLolo$Target.SYBR)

table(lolo$EH_ID, lolo$tissue, lolo$Target.SYBR)

