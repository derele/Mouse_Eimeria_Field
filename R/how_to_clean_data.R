## How to clean input data? To update every year...

## ************* ## NB : correct worms with Jenny!!

## MICE TABLE *****************************

## Jaroslav table genotypes 2013 --> 2016
jarda <- read.csv("../raw_data/HIforEH_May2017.csv", na.strings=c(""," ","NA"))

## Manual names uniformisation:
names(jarda)[names(jarda) %in% c("PIN", "Xmap", "Ymap", "BW", "L", "LCd", "Dissection")] <- 
  c("Mouse_ID", "Longitude", "Latitude", "Body_weight", "Body_length", "Tail_length")

jarda$Mouse_ID <- gsub(pattern = "SK", replacement = "SK_",x = jarda$Mouse_ID)
jarda$Dissection <- paste(paste(jarda$Day, jarda$Month, sep = "/"), jarda$Year, sep = "/") # add dissection date

jarda$Dissection <- as.Date(jarda$Dissection, "%d/%m/%Y") 
jarda$Capture <- jarda$Dissection - jarda$Days.in.lab

# remove useless columns
drops <- names(jarda)[c(2:6, 13, 14, 38, 39, 43:53,59, 66:75)]
jarda <- jarda[ , !(names(jarda) %in% drops)]

# Rename homogeneously
names(jarda)[names(jarda) %in% c("Taenia", "Hymenolepis", "Mastophorus", "Rodentolepis", "Mesocestoides", 
                    "Trichuris", "Heterakis", "Aspiculuris", "Syphacia")] <- 
  c("Cysticercus", "Hymenolepis", "Mastophorus_muris", "Hymenolepis_microstoma", "Mesocestoides",
    "Trichuris_muris", "Heterakis_spumosa", "Aspiculuris_tetraptera", "Syphacia_obvelata" )

micelist <- jarda$Mouse_ID

## **********************************************************
## complement data
diss2014 <- read.csv("../raw_data/HZ14_Mice 31-12-14_dissections.csv", na.strings=c(""," ","NA"))
names(diss2014)[names(diss2014) %in% c("PIN", "X_Map", "Y_Map", "BW", "L", "LCd")] <- 
  c("Mouse_ID", "Longitude", "Latitude", "Body_weight", "Body_length", "Tail_length")
diss2014$Mouse_ID <- paste("SK", diss2014$Mouse_ID, sep = "_")

## rename with homogeneity
names(diss2014)[34:37] <- c("Aspiculuris_tetraptera", "Syphacia_obvelata", "Trichuris_muris", "Cysticercus")

diss2014$Dissection <- as.Date(diss2014$Dissection, "%m/%d/%Y") 
diss2014$Capture <- as.Date(diss2014$Arrival, "%m/%d/%Y")

# remove useless columns
drops <- names(diss2014)[c(2:5, 9:11,16,20:33,39)]
diss2014 <- diss2014[ , !(names(diss2014) %in% drops)]

MiceTable <- merge(jarda, diss2014[!(diss2014$Mouse_ID %in% micelist),], all = TRUE)

# how many new mice? +195
micelist <- MiceTable$Mouse_ID

## **********************************************************
# 2015 - part 1
diss2015.1 <- read.csv("../raw_data/Genotypes_Bav2015.csv", na.strings=c(""," ","NA"))
diss2015.1$Mouse_ID <-  gsub(pattern = "SK", replacement = "SK_",x = diss2015.1$PIN)
names(diss2015.1)[names(diss2015.1) %in% c("Xmap", "Ymap")] <- c("Longitude", "Latitude")
diss2015.1 <- diss2015.1[-c(1,2)]

# 2015 - part 2
diss2015.2 <- read.csv("../raw_data/HZ15_Mice_Parasite.csv", na.strings=c(""," ","NA"))

names(diss2015.2)[names(diss2015.2) %in% c("PIN", "X_Map", "Y_Map", "BW", "L", "LCd", "Aspiculuris.tetraptera",
                                       "Syphacia.obvelata", "Trichuris.muris","Taenia.taeniformis")] <- 
  c("Mouse_ID", "Longitude", "Latitude", "Body_weight", "Body_length", "Tail_length",
    "Aspiculuris_tetraptera", "Syphacia_obvelata", "Trichuris_muris", "Cysticercus")
diss2015.2$Mouse_ID <- paste("SK", diss2015.2$Mouse_ID, sep = "_")

diss2015.2$Dissection <- as.Date(diss2015.2$Dissection, "%m/%d/%Y") 
diss2015.2$Capture <- diss2015.2$Dissection - 1

diss2015.2 <- diss2015.2[-c(2:5, 10:12, 16, 18, 19, 23:36,42)]

# 2015 merged
diss2015 <- merge(diss2015.1, diss2015.2, by = "Mouse_ID", all = TRUE)

# clean
diss2015$Sex <- diss2015$Sex.x
diss2015$Sex[is.na(diss2015$Sex)] <- diss2015$Sex.y[is.na(diss2015$Sex)] 

diss2015$Year <- 2015

diss2015$Longitude <- diss2015$Longitude.x
diss2015$Longitude[is.na(diss2015$Longitude)] <- diss2015$Longitude.y[is.na(diss2015$Longitude)]

diss2015$Latitude <- diss2015$Latitude.x
diss2015$Latitude[is.na(diss2015$Latitude)] <- diss2015$Latitude.y[is.na(diss2015$Latitude)]

diss2015 <- diss2015[-c(2:5, 20, 24,25,26, 28, 29)]

MiceTable <- merge(MiceTable, diss2015[!(diss2015$Mouse_ID %in% micelist),], all = TRUE)

MiceTable <- MiceTable[MiceTable$Mouse_ID != "SK_NA",]

micelist <- MiceTable$Mouse_ID

## 2016
diss2016 <- read.csv("../raw_data/HZ16_Mice_18-07-16_dissections.csv", na.strings=c(""," ","NA"))[-c(1:2),]
names(diss2016)[names(diss2016) %in% "ID_mouse"] <- "Mouse_ID"
diss2016$Mouse_ID <- as.character(diss2016$Mouse_ID)
## Add worms 
worms16 <- read.csv("../raw_data/HZ16_Worms.csv", na.strings=c(""," ","NA"))[-11]
## rename with homogeneity
names(worms16) <- c("Mouse_ID", "Syphacia_obvelata", "Aspiculuris_tetraptera", "Mix_Syphacia_Aspiculuris",
                    "Heterakis_spumosa", "Mastophorus_muris", "Trichuris_muris", 
                    "Hymenolepis_microstoma", "Catenotaenia_pusilla", "Cysticercus")
## merge worms and dissection table 2016
diss2016 <- merge(diss2016, worms16, all = TRUE)
# rename for homogeneity
names(diss2016)[names(diss2016) %in% "Ectoparasites"] <- "Flea"

# remove useless
diss2016 <- diss2016[-c(3:5, 11, 16:22)]

diss2016$Capture <- as.Date(diss2016$Capture, "%d.%m.%Y") 
diss2016$Dissection <- as.Date(diss2016$Dissection, "%d.%m.%Y") 

# For the mice of 2016 in micelist, change the worms column in MiceTable (corrected by Jenny)

wormscol1 <- c("Syphacia_obvelata", "Aspiculuris_tetraptera", "Mix_Syphacia_Aspiculuris",
               "Heterakis_spumosa", "Mastophorus_muris", "Trichuris_muris", 
               "Hymenolepis_microstoma", "Catenotaenia_pusilla", "Cysticercus")

wormscol2 <- c("Syphacia_obvelata", "Aspiculuris_tetraptera", "Heterakis_spumosa", 
               "Mastophorus_muris", "Trichuris_muris", "Hymenolepis_microstoma", 
               "Cysticercus",  "Hymenolepis",  "Mesocestoides",   "Oxyurids")

# extract and correct
extract <- MiceTable[MiceTable$Mouse_ID %in% diss2016$Mouse_ID,]
#remove worms col
extract <- extract[ , !(names(extract) %in% wormscol2)]
# inject correct columns
extract <- merge(extract, diss2016[names(diss2016) %in% c("Mouse_ID",wormscol1)], by = "Mouse_ID", all = TRUE)

# re-inject, corrected
MiceTable <- merge(MiceTable[!(MiceTable$Mouse_ID %in% diss2016$Mouse_ID),], extract, all = TRUE)

## **********************************************************
## 2017
diss2017 <- read.csv(file = "../raw_data/HZ17_September_Mice_Dissection.csv", na.strings=c(""," ","NA"))
names(diss2017)[names(diss2017) == "Ectoparasites"] <- "Flea"      

## Add worms 
worms17 <- read.csv("../raw_data/HZ17_Worms.csv", na.strings=c(""," ","NA"))[-c(13:15)]

## merge worms and dissection table 2017
diss2017 <- merge(diss2017, worms17)

# Remove address and code 
drops <- c("Address","Code",names(diss2017)[17:26]) ; diss2017 <- diss2017[ , !(names(diss2017) %in% drops)]

diss2017$Capture <- as.Date(diss2017$Capture, "%d.%m.%Y") 
diss2017$Dissection <- as.Date(diss2017$Dissection, "%d.%m.%Y") 

# Final merging :) :) :)
MiceTable <- merge(MiceTable, diss2017, all = TRUE)

## Uniformisation:
levels(MiceTable$Sex)[levels(MiceTable$Sex) %in% c("female", "male")] <- c("F", "M")

write.csv(x = MiceTable, file = "../raw_data/MiceTable_2014to2017.csv", row.names = FALSE)

## TRAPS TABLE *****************************

# 2017
traps2017 <- read.csv(file = "../raw_data/HZ17_Mice_Trap.csv")
traps2017 <- traps2017[-1]
traps2017$Year <- 2017
traps2017[traps2017$Number_mus_caught == "NA (“keiner da”)", "Number_mus_caught"] <- 0
traps2017$Number_mus_caught <- as.numeric(as.character(traps2017$Number_mus_caught))

# 2016
traps2016 <- read.csv("../raw_data/Cleaned_HMHZ_2016_All.csv")
traps2016$Year <- 2016
# remove garbage
drops <- c("X","Got.the.sheet.", "location", "area", "Y_N", "MperT") 
traps2016 <- traps2016[ , !(names(traps2016) %in% drops)]

# uniformise names
names(traps2016) <- c("Address", "Latitude", "Longitude", "Date_set", "Time_set", "Number_traps_set", "Date_collect", 
    "Time_collect","Number_mus_caught", "People", "Number_rodents_caught", "Year")
  
# before 2016, only localisations where mice where actually caught
trapsbefore17 <- MiceTable[names(MiceTable) %in% c("Longitude", "Latitude", "Year")]
trapsbefore17 <- trapsbefore17[!trapsbefore17$Year %in% c(2017, 2016), ]
# each line is one mouse in the dissection table
trapsbefore17$Number_mus_caught <- 1

trapsbefore17 <- aggregate(Number_mus_caught ~ Latitude + Longitude + Year, FUN = sum, data = trapsbefore17)

trapsbefore17$People <- NA; trapsbefore17$Address <- NA; trapsbefore17$Date_set <- NA; trapsbefore17$Time_set <-NA
trapsbefore17$Number_traps_set <- NA; trapsbefore17$Date_collect <- NA; trapsbefore17$Time_collect <- NA;
trapsbefore17$Number_rodents_caught <- NA; trapsbefore17$People <- NA

TrapTable_2014to2017 <- rbind(rbind(trapsbefore17, traps2016), traps2017)

write.csv(x = TrapTable_2014to2017, file = "../raw_data/TrapTable_2014to2017.csv", row.names = FALSE)

