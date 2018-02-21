## How to clean input data? To update every year...
source("HMHZ_Functions.R")
library(data.table)
library(ggmap)
## ************* ## NB : correct worms with Jenny!!

## MICE TABLE *****************************
## Jaroslav table genotypes 2013 --> 2017
HIJardaTable <- read.csv("../raw_data/EmanuelData.csv", na.strings=c(""," ","NA"), stringsAsFactors = FALSE)
HIJardaTable <- HIJardaTable[!names(HIJardaTable) == "X"]

## Manual names uniformisation:
setnames(HIJardaTable,
         old = c("PIN", "X_Longit", "Y_Latit"), 
         new = c("Mouse_ID", "Longitude", "Latitude"))

HIJardaTable$Mouse_ID <- gsub(pattern = "SK", replacement = "SK_",x = HIJardaTable$Mouse_ID)

# Let's remove the embryos (sadly, no interest for our parasitic studies...)
HIJardaTable <- HIJardaTable[sapply(HIJardaTable$Mouse_ID, nchar) <= 7,]

# A first map to be happy about the hard work :)
HI.map(HIJardaTable)

# How many samples from Brandenburg do we have the HI for per year?
table(HIJardaTable$Year)

# CHECK DUPLICATE
HIJardaTable$Mouse_ID[duplicated(HIJardaTable$Mouse_ID)]

# Complement data with extra HI calculated 
# 2014 : no extra HI; 2015 : extra HIs :)

## 2015 part 1
diss2015.1 <- read.csv("../raw_data/Genotypes_Bav2015.csv", na.strings=c(""," ","NA"), stringsAsFactors = FALSE)

# Manual names uniformisation
setnames(diss2015.1,
         old = c("PIN", "Xmap", "Ymap"), 
         new = c("Mouse_ID", "Longitude", "Latitude"))
diss2015.1$Mouse_ID <-  gsub(pattern = "SK", replacement = "SK_", x = diss2015.1$Mouse_ID)

# Add extra mice to full miceTable
extraMice <- diss2015.1$Mouse_ID[!diss2015.1$Mouse_ID %in% HIJardaTable$Mouse_ID]
miceTable <- merge(HIJardaTable, diss2015.1[diss2015.1$Mouse_ID %in% extraMice,], all = T)
miceTable$Mouse_ID[duplicated(miceTable$Mouse_ID)]

## 2015 part 2
diss2015.2 <- read.csv("../raw_data/HZ15_Mice_Parasite.csv", na.strings=c(""," ","NA"), stringsAsFactors = FALSE)

# remove transported mice
diss2015.2 <- diss2015.2[!is.na(diss2015.2$PIN),]

# Manual names uniformisation
setnames(diss2015.2,
         old = c("X", "PIN", "X_Map", "Y_Map", "BW", "L", "LCd", 
                 "Code", "State", "Locality", "Sex",
                 "Aspiculuris.tetraptera", "Syphacia.obvelata", "Trichuris.muris","Taenia.taeniformis"),
         new = c("Notes", "Mouse_ID", "Longitude.extra", "Latitude.extra", "Body_weight", "Body_length", "Tail_length", 
                 "Code2015", "State2015", "Locality2015", "Sex2015",
                 "Aspiculuris_tetraptera", "Syphacia_obvelata", "Trichuris_muris", "Cysticercus"))
diss2015.2$Mouse_ID <- paste0(diss2015.2$ID, "_", diss2015.2$Mouse_ID)

# Fill up transect info from HZ_BAV in miceTable
miceWithoutTransect <- miceTable$Mouse_ID[is.na(miceTable$Transect)]
miceWithTransectHZ_BAV <- diss2015.2$Mouse_ID[diss2015.2$Mouse_ID %in% miceWithoutTransect]
miceTable$Transect[miceTable$Mouse_ID %in% miceWithTransectHZ_BAV] <- "HZ_BAV"

# some have their transect lost... TODO LATER
miceTable$Mouse_ID[is.na(miceTable$Transect)]

# Add extra dissection info to full miceTable
extraInfoDissection <- diss2015.2[diss2015.2$Mouse_ID %in% miceTable$Mouse_ID,]

miceTable <- merge(miceTable, extraInfoDissection, all = T)

# Fill if Lat/lon missing, from 2015 added DF
miceTable$Latitude[is.na(miceTable$Latitude)] <- miceTable$Latitude.extra[is.na(miceTable$Latitude)]
miceTable$Longitude[is.na(miceTable$Longitude)] <- miceTable$Longitude.extra[is.na(miceTable$Longitude)]

###################### 2016
diss2016 <- read.csv("../raw_data/HZ16_Mice_18-07-16_dissections.csv", na.strings=c(""," ","NA"))[-c(1:2),]
names(diss2016)[names(diss2016) %in% "ID_mouse"] <- "Mouse_ID"
diss2016$Mouse_ID <- as.character(diss2016$Mouse_ID)

## Add worms 
worms16 <- read.csv("../raw_data/HZ16_Worms.csv", na.strings=c(""," ","NA"))[-11]
## rename with homogeneity
names(worms16) <- c("Mouse_ID", "Syphacia_obvelata", "Aspiculuris_tetraptera", "Mix_Syphacia_Aspiculuris",
                    "Heterakis_spumosa", "Mastophorus_muris", "Trichuris_muris", 
                    "Hymenolepis_microstoma", "Catenotaenia_pusilla", "Cysticercus")

worms16$Oxyurids <- rowSums(worms16[c("Aspiculuris_tetraptera",
                                      "Syphacia_obvelata", 
                                      "Mix_Syphacia_Aspiculuris")])

## merge worms and dissection table 2016
diss2016 <- merge(diss2016, worms16, all = TRUE)

# rename for homogeneity
names(diss2016)[names(diss2016) %in% "Ectoparasites"] <- "Flea"

diss2016$Capture <- as.Date(diss2016$Capture, "%d.%m.%Y") 
diss2016$Dissection <- as.Date(diss2016$Dissection, "%d.%m.%Y") 

# merge
NewmiceTable <- merge(miceTable, diss2016, 
                      by = c("Mouse_ID"), all = T)

NewmiceTable <- fillGapsAfterMerge(NewmiceTable)

## **********************************************************
## 2017
diss2017 <- read.csv(file = "../raw_data/HZ17_September_Mice_Dissection.csv", na.strings=c(""," ","NA"), stringsAsFactors = F)

# correction excel bullshit
diss2017$Feces_weight <- as.numeric(as.character(diss2017$Feces_weight))

diss2017$Feces_weight[diss2017$Feces_weight > 100 & !is.na(diss2017$Feces_weight)] <-
  diss2017$Feces_weight[diss2017$Feces_weight > 100 & !is.na(diss2017$Feces_weight)] / 1000

names(diss2017)[names(diss2017) == "Ectoparasites"] <- "Flea"      

## Add worms 
worms17 <- read.csv2("../raw_data/HZ17_September_Mice_Dissection_Jen_final.csv")

names(worms17)[names(worms17) %in% "Mesocestoides"] <- "Taenia_martis"

## merge worms and dissection table 2017
diss2017 <- merge(diss2017, worms17[names(worms17) %in% c("Mouse_ID", "Hymenolepis_microstoma", "Hymenolepis_diminiuta",
                                                          "Catenotaenia_pusilla", "Cysticercus", "Taenia_martis", "Mastophorus_muris",
                                                          "Trichuris_muris","Heterakis_spumosa", "Syphacia_obvelata", "Aspiculuris_tetraptera",
                                                          "Heligmosomoides_polygurus")])

# Remove useless columns
drops <- c("State", "Address","Code", "Dissectors", "Spleen_mass", "Left_Testis_mass" ,"Right_Testis_mass", "Left.epididymis.weight",
           "Seminal.vesicle.weight", "Left.ovarium.weight", "Right.ovarium.weight", "Total", "Embryo_left", "Embryo_right")   

diss2017 <- diss2017[ , !(names(diss2017) %in% drops)]

diss2017$Oxyurids <- rowSums(diss2017[c("Syphacia_obvelata", "Aspiculuris_tetraptera")])

# check
which(duplicated(diss2017$Mouse_ID))

# merge
MiceTable <- merge(NewmiceTable, diss2017, 
                      by = c("Mouse_ID"), all = T)

MiceTable <- fillGapsAfterMerge(MiceTable)

## Uniformisation:
levels(MiceTable$Sex)[levels(MiceTable$Sex) %in% c("female", "male")] <- c("F", "M")

# Proper names
names(MiceTable)[names(MiceTable) %in% "Oxyurids"] <- "Oxyuridea"
names(MiceTable)[names(MiceTable) %in% "Mesocestoides"] <- "Mesocestoides_sp."
names(MiceTable)[names(MiceTable) %in% "Cysticercus"] <- "Taenia_taeniformis"

# Remove if no HI
MiceTable <- MiceTable[!is.na(MiceTable$HI),]

# Final plot!
HI.map(MiceTable)











############ CLEAN
## Groups
MiceTable$Hymenolepis[is.na(MiceTable$Hymenolepis)] <- 
  MiceTable$Hymenolepis_microstoma[is.na(MiceTable$Hymenolepis)] +
  MiceTable$Hymenolepis_diminiuta[is.na(MiceTable$Hymenolepis)] +
  MiceTable$Hymenolepis_diminiuta[is.na(MiceTable$Hymenolepis)]

MiceTable$Tapeworms <- 
  rowSums(MiceTable[c("Taenia_martis", "Taenia_taeniformis",
                    "Mesocestoides_sp.","Catenotaenia_pusilla",
                    "Hymenolepis")], na.rm = T)

## Big groups WATWM = Tapeworms, Whipworms (t.muris), oxyuridae, m.muris
WormsDF <- MiceTable[c("Mouse_ID", "Year", "Tapeworms",
                       "Trichuris_muris", "Oxyuridea", "Mastophorus_muris")]

library(reshape)
WormsDF <- melt(WormsDF, id = c("Mouse_ID", "Year"))

WormsDF$value <- as.numeric(as.character(WormsDF$value))

ggplot(data=WormsDF, aes(x = variable, y=log10(value))) +
  geom_violin(aes(fill = variable))  +
  geom_jitter(size = 0.5, alpha = .8) +
  theme_classic() +
  facet_wrap( ~ Year, nrow = 2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position="none")

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

