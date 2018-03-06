## How to clean input data? To update every year...
source("HMHZ_Functions.R")
library(data.table)
library(ggmap)
library(reshape)
## ************* ## NB : correct worms with Jenny!!

## MICE TABLE *****************************
## Jaroslav table genotypes 2014 --> 2017
HIJardaTable <- read.csv("../raw_data/EmanuelData.csv", 
                         na.strings=c(""," ","NA"), stringsAsFactors = FALSE)
HIJardaTable <- HIJardaTable[!names(HIJardaTable) == "X"]
HIJardaTable <- HIJardaTable[!HIJardaTable$Year %in% c(2010, 2011),]

## Manual names uniformisation:
setnames(HIJardaTable,
         old = c("PIN", "X_Longit", "Y_Latit"), 
         new = c("Mouse_ID", "Longitude", "Latitude"))

HIJardaTable$Mouse_ID <- gsub(pattern = "SK", replacement = "SK_",x = HIJardaTable$Mouse_ID)

# Let's remove the embryos (sadly, no interest for our parasitic studies...)
HIJardaTable <- HIJardaTable[sapply(HIJardaTable$Mouse_ID, nchar) <= 7,]

# How many samples from Brandenburg do we have the HI for per year?
table(HIJardaTable$Year)

# CHECK DUPLICATE
HIJardaTable$Mouse_ID[duplicated(HIJardaTable$Mouse_ID)]

# Complement data with previous tables (2014, 2015) 
diss2014.1 <- read.csv("../raw_data/HZ14_Mice 31-12-14_dissections.csv",
                     na.strings=c(""," ","NA"), stringsAsFactors = FALSE)

# Homogenize : Mouse_ID, lat, lon, worms
diss2014.1$Mouse_ID <- paste0(diss2014.1$ID, "_", diss2014.1$PIN)
setnames(diss2014.1,
         old = c("X_Map", "Y_Map", 
                 "Aspiculuris.tetraptera", "Syphacia.obvelata",
                 "Trichuris.muris", "Taenia.taeniformis"), 
         new = c("Longitude", "Latitude",
                 "Aspiculuris_tetraptera", "Syphacia_obvelata",
                 "Trichuris_muris", "Taenia_taeniformis"))

# Genotypes 2014
diss2014.2 <- read.csv("../raw_data/HZ14_Mice 31-12-14_genotypes.csv")
diss2014.2$Mouse_ID <- paste0(diss2014.2$ID, "_", diss2014.2$PIN)

# Calculta HIs
diss2014.2 <- get.HI.full(diss2014.2)
setnames(diss2014.2, old = "HI.calculated", new = "HI")

# Merge & complete
diss2014 <- merge(diss2014.1, diss2014.2,  
                         by = c("Mouse_ID"), all = T)
diss2014 <- fillGapsAfterMerge(diss2014)

# Merge & complete
mergedMiceTable <- merge(HIJardaTable, diss2014,  
                         by = c("Mouse_ID"), all = T)
mergedMiceTable <- fillGapsAfterMerge(mergedMiceTable)

## 2015 part 1
diss2015.1 <- read.csv("../raw_data/Genotypes_Bav2015.csv", 
                       na.strings=c(""," ","NA"), stringsAsFactors = FALSE)

# Manual names uniformisation
setnames(diss2015.1,
         old = c("PIN", "Xmap", "Ymap"), 
         new = c("Mouse_ID", "Longitude", "Latitude"))
diss2015.1$Mouse_ID <-  gsub(pattern = "SK", replacement = "SK_", x = diss2015.1$Mouse_ID)

# Merge & complete
mergedMiceTable <- merge(mergedMiceTable, diss2015.1,  
                         by = c("Mouse_ID"), all = T)
mergedMiceTable <- fillGapsAfterMerge(mergedMiceTable)

## 2015 part 2
diss2015.2 <- read.csv("../raw_data/HZ15_Mice_Parasite.csv", 
                       na.strings=c(""," ","NA"), stringsAsFactors = FALSE)

# remove transported mice
diss2015.2 <- diss2015.2[!is.na(diss2015.2$PIN),]

# Manual names uniformisation
setnames(diss2015.2,
         old = c("X", "PIN", "X_Map", "Y_Map", "BW", "L", "LCd", 
                 "Aspiculuris.tetraptera", "Syphacia.obvelata", 
                 "Trichuris.muris","Taenia.taeniformis"),
         new = c("Notes", "Mouse_ID", "Longitude.extra", "Latitude.extra", "Body_weight", "Body_length", "Tail_length", 
                 "Aspiculuris_tetraptera", "Syphacia_obvelata", 
                 "Trichuris_muris", "Taenia_taeniformis"))
diss2015.2$Mouse_ID <- paste0(diss2015.2$ID, "_", diss2015.2$Mouse_ID)

# Merge & complete
mergedMiceTable <- merge(mergedMiceTable, diss2015.2,  
                         by = c("Mouse_ID"), all = T)
mergedMiceTable <- fillGapsAfterMerge(mergedMiceTable)

# some have their transect lost... TODO LATER
mergedMiceTable$Mouse_ID[is.na(mergedMiceTable$Transect)]

###################### 2016
diss2016 <- read.csv("../raw_data/HZ16_Mice_18-07-16_dissections.csv", 
                     na.strings=c(""," ","NA"), stringsAsFactors = F)[-c(1:2),]
names(diss2016)[names(diss2016) %in% "ID_mouse"] <- "Mouse_ID"
diss2016$Mouse_ID <- as.character(diss2016$Mouse_ID)

## Add worms 
worms16 <- read.csv("../raw_data/HZ16_Worms.csv", 
                    na.strings=c(""," ","NA"), stringsAsFactors = F)[-11]
## rename with homogeneity
names(worms16) <- c("Mouse_ID", "Syphacia_obvelata", "Aspiculuris_tetraptera", "Mix_Syphacia_Aspiculuris",
                    "Heterakis_spumosa", "Mastophorus_muris", "Trichuris_muris", 
                    "Hymenolepis_microstoma", "Catenotaenia_pusilla", "Cysticercus")

## merge worms and dissection table 2016
diss2016 <- merge(diss2016, worms16, all = TRUE)

# rename for homogeneity
names(diss2016)[names(diss2016) %in% "Ectoparasites"] <- "Flea"
diss2016$Capture <- as.Date(diss2016$Capture, "%d.%m.%Y") 
diss2016$Dissection <- as.Date(diss2016$Dissection, "%d.%m.%Y") 

# merge
mergedMiceTable <- merge(mergedMiceTable, diss2016, 
                      by = c("Mouse_ID"), all = T)

mergedMiceTable <- fillGapsAfterMerge(mergedMiceTable)

## **********************************************************
## 2017
diss2017 <- read.csv(file = "../raw_data/HZ17_September_Mice_Dissection.csv",
                     na.strings=c(""," ","NA"), stringsAsFactors = F)

# correction excel bullshit
diss2017$Feces_weight <- as.numeric(as.character(diss2017$Feces_weight))

diss2017$Feces_weight[diss2017$Feces_weight > 100 & !is.na(diss2017$Feces_weight)] <-
  diss2017$Feces_weight[diss2017$Feces_weight > 100 & !is.na(diss2017$Feces_weight)] / 1000

names(diss2017)[names(diss2017) == "Ectoparasites"] <- "Flea"      

## Add worms 
worms17 <- read.csv2("../raw_data/HZ17_September_Mice_Dissection_Jen_final.csv",
                     stringsAsFactors = F)

names(worms17)[names(worms17) %in% "Mesocestoides"] <- "Taenia_martis"

## merge worms and dissection table 2016
diss2017 <- merge(diss2017, worms17, 
                  by = "Mouse_ID", all = TRUE)
diss2017 <- fillGapsAfterMerge(diss2017)

# merge
mergedMiceTable <- merge(mergedMiceTable, diss2017, 
                      by = c("Mouse_ID"), all = T)

mergedMiceTable <- fillGapsAfterMerge(mergedMiceTable)

## Uniformisation
mergedMiceTable$Sex[grep("female*.", mergedMiceTable$Sex)] <- "F"
mergedMiceTable$Sex[grep("male*.", mergedMiceTable$Sex)] <- "M"

# Remove if no HI
mergedMiceTable <- mergedMiceTable[!is.na(mergedMiceTable$HI),]

############ Worms ############
## in WATWM dataset : Hymenolepis, Taenia, Rodentolepis, Mesocestoides,
## Calodium, Mastophorus, Trichuris, Heterakis, Aspiculuris+Syphacia

# Hymenolepis
mergedMiceTable$Hymenolepis <- rowSums(
  mergedMiceTable[c("Hymenolepis_microstoma", "Hymenolepis_diminiuta")], 
  na.rm = T)
mergedMiceTable$Hymenolepis[with(mergedMiceTable, 
                                 is.na(mergedMiceTable["Hymenolepis_microstoma"]) &  
                                   is.na(mergedMiceTable["Hymenolepis_diminiuta"]))] <- NA

# Taenia
mergedMiceTable$Taenia <- rowSums(
  mergedMiceTable[c("Taenia_martis", "Taenia_taeniformis", 
                    "Catenotaenia_pusilla", "Cysticercus")], 
  na.rm = T)
mergedMiceTable$Hymenolepis[with(mergedMiceTable, 
                                 is.na(mergedMiceTable["Taenia_martis"]) &  
                                   is.na(mergedMiceTable["Taenia_taeniformis"]) &
                                   is.na(mergedMiceTable["Catenotaenia_pusilla"]) &
                                   is.na(mergedMiceTable["Cysticercus"]))] <- NA

# Aspiculuris_Syphacia
mergedMiceTable$Aspiculuris_Syphacia <- rowSums(
  mergedMiceTable[c("Syphacia_obvelata", "Aspiculuris_tetraptera", "Mix_Syphacia_Aspiculuris")], 
  na.rm = T)
mergedMiceTable$Aspiculuris_Syphacia[with(mergedMiceTable, 
                                 is.na(mergedMiceTable["Syphacia_obvelata"]) &  
                                   is.na(mergedMiceTable["Aspiculuris_tetraptera"]) &
                                   is.na(mergedMiceTable["Mix_Syphacia_Aspiculuris"]))] <- NA

# Trichuris
mergedMiceTable$Trichuris <- mergedMiceTable$Trichuris_muris

# Heterakis
mergedMiceTable$Heterakis <- mergedMiceTable$Heterakis_spumosa

# Mastophorus
mergedMiceTable$Mastophorus <- mergedMiceTable$Mastophorus_muris

## Dataframe to plot worms
WormsDF <- mergedMiceTable[c("Mouse_ID", "Year", 
                             "Hymenolepis", "Taenia", "Aspiculuris_Syphacia", 
                             "Trichuris", "Heterakis", "Mastophorus")]
WormsDF <- melt(WormsDF, id = c("Mouse_ID", "Year"))

WormsDF$value <- as.numeric(as.character(WormsDF$value))

ggplot(data=WormsDF, aes(x = variable, y=log10(value))) +
  geom_violin(aes(fill = variable))  +
  geom_jitter(size = 0.5, width = .2, alpha = .8) +
  theme_classic() +
  facet_wrap( ~ Year, nrow = 2) +
  theme(text = element_text(size = 15),
        axis.text = element_text(angle = 45, hjust = 1))+
  theme(legend.position="none")
## TODO 2014 and 2015 worms, not correct!!

# Compare to WATWM data!
# WATWM <- read.csv("https://raw.githubusercontent.com/alicebalard/Parasite_Load/refactor/examples/Reproduction_WATWM/EvolutionFinalData.csv")
# WormsWATWMDF <- melt(WATWM, id = c("Individual_mouse", "Latitude", "Longtitude",
#                                    "Sex", "HI"))
# 
# WormsWATWMDF$value <- as.numeric(as.character(WormsWATWMDF$value))
# 
# ggplot(data = WormsWATWMDF, aes(x = variable, y=log10(value))) +
#   geom_violin(aes(fill = variable))  +
#   geom_jitter(size = 0.5, width = .2, alpha = .8) +
#   theme_classic() +
#   theme(text = element_text(size = 15),
#         axis.text = element_text(angle = 45, hjust = 1))+
#   theme(legend.position="none")

# Final cleaning, and save!
mergedMiceTable$Longitude <- as.numeric(mergedMiceTable$Longitude)
mergedMiceTable$Latitude <- as.numeric(mergedMiceTable$Latitude)

# Remove far East samples
mergedMiceTable <- mergedMiceTable[mergedMiceTable$Longitude < 18,]

# How many samples from Brandenburg do we have the HI for per year?
table(mergedMiceTable$Year, mergedMiceTable$Transect)

# correct body length/weight
mergedMiceTable$Body_length <- as.numeric(mergedMiceTable$Body_length)

mergedMiceTable$Body_weight <- as.numeric(mergedMiceTable$Body_weight)
# Manual correction
mergedMiceTable$Body_weight[mergedMiceTable$Body_weight >= 200 & !is.na(mergedMiceTable$Body_weight)]  <- 
mergedMiceTable$Body_weight[mergedMiceTable$Body_weight >= 200 & !is.na(mergedMiceTable$Body_weight)] / 1000

# Body condition index as log body mass/log body length (Hayes et al. 2014)
miceTable$BCI <- log(miceTable$Body_weight) / log(miceTable$Body_length)

########## Write out ##########
write.csv(x = mergedMiceTable, file = "../raw_data/MiceTable_2014to2017.csv", row.names = FALSE)

# ## plot for following over the years
# HI.map(df = mergedMiceTable[mergedMiceTable$Year == 2017,], margin = .2)
# HI.map(df = mergedMiceTable, margin = .2)
# 
# df <- data.frame(HI = mergedMiceTable$HI, year = as.factor(mergedMiceTable$Year))
# 
# ggplot(data = df, 
#        aes(x = year, y = HI, color = HI))+
#   geom_violin() +
#   geom_jitter() +
#   scale_color_gradient("Hybrid\nindex", high="red",low="blue") +
#   theme_bw()
# 
# ## Distribution of HI just for BR transect
# dfBR <- mergedMiceTable[mergedMiceTable$Transect == "HZ_BR",]
# df <- data.frame(HI = dfBR$HI, year = as.factor(dfBR$Year))
# 
# ggplot(data = df, 
#        aes(x = year, y = HI, color = HI))+
#   geom_violin() +
#   geom_jitter() +
#   scale_color_gradient("Hybrid\nindex", high="red",low="blue") +
#   theme_bw()
# 
# # Worms distribution over HI
# mergedMiceTable$ALLWORMS <- rowSums(
#   mergedMiceTable[c("Mastophorus", "Heterakis", "Trichuris", "Aspiculuris_Syphacia",
#                     "Taenia", "Hymenolepis")], 
#   na.rm = T)
# 
# ggplot(data = mergedMiceTable, 
#        aes(x = HI, y = log10(ALLWORMS), color = HI))+
#   geom_jitter(aes(shape = as.factor(Year)), size = 3) +
#   geom_smooth() +
#   theme(text = element_text(size = 15)) +
#   scale_color_gradient("Hybrid\nindex", high="red",low="blue") +
#   theme_bw()
# 
# ggplot(data = mergedMiceTable, 
#        aes(x = HI, y = log10(Trichuris), color = HI))+
#   geom_jitter(aes(shape = as.factor(Year)), size = 3) +
#   geom_smooth() +
#   theme(text = element_text(size = 15)) +
#   scale_color_gradient("Hybrid\nindex", high="red",low="blue") +
#   theme_bw()
# 
# ggplot(data = mergedMiceTable, 
#        aes(x = HI, y = log10(Aspiculuris_Syphacia), color = HI))+
#   geom_jitter(aes(shape = as.factor(Year)), size = 3) +
#   geom_smooth() +
#   theme(text = element_text(size = 20)) +
#   scale_color_gradient("Hybrid\nindex", high="red",low="blue") +
#   theme_bw()
# 
# Lorenzo <- read.csv("../raw_data/Eimeria_detection/Eimeria_oocysts_2017_Lorenzo.csv")
# LorenzoFull <- merge(Lorenzo, mergedMiceTable, by = "Mouse_ID")
# 
# ggplot(data = LorenzoFull, 
#        aes(x = HI, y = log10(OPG), color = HI))+
#   geom_jitter() +
#   geom_smooth() +
#   theme(text = element_text(size = 20)) +
#   scale_color_gradient("Hybrid\nindex", high="red",low="blue") +
#   theme_bw()
# 
