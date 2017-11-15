## How to clean input data? To update every year...

## ************* ## NB : correct worms with Jenny!!

## MICE TABLE *****************************

## Jaroslav table genotypes 2013 --> 2016
jarda <- read.csv("../raw_data/HIforEH_May2017.csv")

## Manual names uniformisation:
names(jarda)[names(jarda) %in% c("PIN", "Xmap", "Ymap", "BW", "L", "LCd")] <- 
  c("Mouse_ID", "Longitude", "Latitude", "Body_weight", "Body_length", "Tail_length")

jarda$Mouse_ID <- gsub(pattern = "SK", replacement = "SK_",x = jarda$Mouse_ID)
jarda$Dissection <- paste(paste(jarda$Day, jarda$Month, sep = "/"), jarda$Year, sep = "/") # add dissection date

## **********************************************************
## Which are the "forgotten" samples, with no mouse genotype?
diss2014 <- read.csv("../raw_data/HZ14_Mice 31-12-14_dissections.csv")
names(diss2014)[names(diss2014) %in% c("PIN", "X_Map", "Y_Map", "BW", "L", "LCd")] <- 
  c("Mouse_ID", "Longitude", "Latitude", "Body_weight", "Body_length", "Tail_length")
diss2014$Mouse_ID <- paste("SK", diss2014$Mouse_ID, sep = "_")

## 2015 with genotypes
diss2015gen <- read.csv("../raw_data/Genotypes_Bav2015.csv")
names(diss2015gen)[names(diss2015gen) %in% c("PIN", "Xmap", "Ymap")] <- 
  c("Mouse_ID", "Longitude", "Latitude")
diss2015gen$Mouse_ID <-  gsub(pattern = "SK", replacement = "SK_",x = diss2015gen$Mouse_ID)
diss2015gen$Body_weight <- NA; diss2015gen$Body_length <- NA; diss2015gen$Tail_length <- NA

## all 2015 dissections
diss2015 <- read.csv("../raw_data/HZ15_Mice_Parasite.csv")
names(diss2015)[names(diss2015) %in% c("PIN", "X_Map", "Y_Map", "BW", "L", "LCd")] <- 
  c("Mouse_ID", "Longitude", "Latitude", "Body_weight", "Body_length", "Tail_length")
diss2015$Mouse_ID <- paste("SK", diss2015$Mouse_ID, sep = "_")

## 2015 without genotypes
diss2015nogen <- diss2015[diss2015$Mouse_ID %in% diss2015gen$Mouse_ID == FALSE,]

## full 2015
diss2015 <- merge(diss2015gen, diss2015nogen, all = TRUE)

## 2016
diss2016 <- read.csv("../raw_data/HZ16_Mice_18-07-16_dissections.csv")[-c(1:2),]
names(diss2016)[names(diss2016) %in% "ID_mouse"] <- "Mouse_ID"
diss2016$Mouse_ID <- as.character(diss2016$Mouse_ID)

## All years
MiceTable_14to16 <- merge(merge(diss2014, diss2015, all = TRUE), diss2016, all = TRUE)

## **********************************************************
## Begin from jarda's table as a base, then complete with missing present in MiceTable_14to16
MiceTable_14to16 <- merge(jarda, MiceTable_14to16[MiceTable_14to16$Mouse_ID %in% jarda$Mouse_ID == FALSE,], all = TRUE)

## **********************************************************
## Add 2017
Mice_17 <- read.csv(file = "../raw_data/HZ17_September_Mice_Dissection.csv")
names(Mice_17)[names(Mice_17) == "Ectoparasites"] <- "Flea"      

# Remove address and code 
drops <- c("Address","Code") ; Mice_17 <- Mice_17[ , !(names(Mice_17) %in% drops)]

# Merge 2017 and MiceTable_14to16
MiceTable_14to17 <- merge(MiceTable_14to16, Mice_17, all = TRUE)

# remove mistakes
MiceTable_14to17 <- MiceTable_14to17[MiceTable_14to17$Mouse_ID != "SK_NA", ]

# Check for unicity:
which(table(MiceTable_14to17$Mouse_ID) != 1)
## **********************************************************

## We keep the following:
#MiceTable_14to17 <- MiceTable_14to17[c("Transect", "State", "PIN", "Sex", "Xmap", "Ymap", "Year", "Month", "Day", "mtBamH", "Zfy2", "SRY1", "Y", "X332", "X347", "X65",
 #                        "Tsx", "Btk", "Syap1", "Es1C", "Gpd1C", "Idh1C", "MpiC", "NpC", "Sod1C", "Days.in.lab", "BW", "L", "LCd", "Taenia",
  #                       "Hymenolepis", "Mastophorus", "Rodentolepis", "Mesocestoides", "Trichuris", "Heterakis", "Aspiculuris",
   #                      "Syphacia", "Oxyurids", "Flea", "HI_NLoci", "HI")]

# remove useless factors to final df
#drops <- c("Dissectors", "Spleen_mass", "Left_Testis_mass", "Right_Testis_mass","Left.epididymis.weight", "Seminal.vesicle.weight", "Left.ovarium.weight", "Right.ovarium.weight", "Total",
 #          "Embryo_left", "Embryo_right"); Mice_14to17 <- Mice_14to17[ , !(names(Mice_14to17) %in% drops)]

## Uniformisation:
levels(MiceTable_14to17$Sex)[levels(MiceTable_14to17$Sex) %in% c("female", "male")] <- c("F", "M")

write.csv(x = MiceTable_14to17, file = "../raw_data/MiceTable_2014to2017.csv", row.names = FALSE)

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
drops <- c("X","Got.the.sheet.", "location", "area", "Y_N", "MperT") ; traps2016 <- traps2016[ , !(names(traps2016) %in% drops)]

# uniformise names
names(traps2016) <- c("Address", "Longitude", "Latitude",  "Date_set", "Time_set", "Number_traps_set", "Date_collect", 
    "Time_collect","Number_mus_caught", "People", "Number_rodents_caught", "Year")
  
# before 2016, only localisations where mice where actually caught
trapsbefore17 <- MiceTable_14to17[names(MiceTable_14to17) %in% c("Latitude", "Longitude", "Year")]
trapsbefore17 <- trapsbefore17[!trapsbefore17$Year %in% c(2017, 2016), ]
# each line is one mouse in the dissection table
trapsbefore17$Number_mus_caught <- 1

trapsbefore17 <- aggregate(Number_mus_caught ~ Latitude + Longitude + Year, FUN = sum, data = trapsbefore17)

trapsbefore17$People <- NA; trapsbefore17$Address <- NA; trapsbefore17$Date_set <- NA; trapsbefore17$Time_set <-NA
trapsbefore17$Number_traps_set <- NA; trapsbefore17$Date_collect <- NA; trapsbefore17$Time_collect <- NA;
trapsbefore17$Number_rodents_caught <- NA; trapsbefore17$People <- NA

TrapTable_2014to2017 <- rbind(rbind(trapsbefore17, traps2016), traps2017)

write.csv(x = TrapTable_2014to2017, file = "../raw_data/TrapTable_2014to2017.csv", row.names = FALSE)