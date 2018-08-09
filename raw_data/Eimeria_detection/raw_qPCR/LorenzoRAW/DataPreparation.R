lolo <- read.csv("Schreibtisch/Students/Lorenzo/Lorenzo/CSVFiles/TotalLorenzo.csv", 
                 na.strings = c("NA", "", " ", "-"))

# Column file name
names(lolo)[names(lolo) == "QPCR01.06.2018.XLS.csv"] <- "fileName"

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

# Add tissue and EH_ID
x <-strsplit(as.character(lolo$Name), "_", 1)

lolo$tissue <- sapply( x, "[", 1)

lolo$EH_ID <- paste0("AA_", sapply( x, "[", 3))

table(lolo$fileName) #QPCR03.05.2018.XLS.csv not complete!!!

# Count triplicates
library(plyr)
summary <- plyr::count(lolo[, c("EH_ID", "tissue", "Target.SYBR", "fileName")])

toCheck <- plyr::count(lolo[, c("EH_ID", "tissue", "Target.SYBR")])

length(which(toCheck$freq != 3))




table(lolo$EH_ID, lolo$tissue, lolo$Target.SYBR))

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

