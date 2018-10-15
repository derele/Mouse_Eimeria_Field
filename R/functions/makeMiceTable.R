makeMiceTable <- function(pathToMyData = "data/Field_data/",
                          pathtoHMHZfunction= "R/functions/HMHZ_Functions.R"){
  
  pathToMyData = "data/Field_data/"
  pathtoHMHZfunction= "R/functions/HMHZ_Functions.R"

  #### NB: all data are locally kept in a folder Data_important
  #### TO BE DEFINED HERE
  
  ## How to clean input data? To update every year...
  source(pathtoHMHZfunction)
  library(data.table)
  library(ggmap)
  library(reshape)
  ## ************* ## NB : correct worms with Jenny!!
  
  ## MICE TABLE *****************************
  ## Jaroslav table genotypes 2014 --> 2017 BEWARE IT LACKS SAMPLES !!!!!
  HIJardaTable <- read.csv(paste0(pathToMyData, "EmanuelData.csv"), 
                           na.strings=c(""," ","NA"), stringsAsFactors = FALSE)
  HIJardaTable <- HIJardaTable[!names(HIJardaTable) == "X"]
  HIJardaTable <- HIJardaTable[!HIJardaTable$Year %in% c(2010, 2011),]
  
  # Check if all rows are NA and delete these rows
  which(!rowSums(!is.na(HIJardaTable)))  
  
  ## Manual names uniformisation:
  setnames(HIJardaTable,
           old = c("PIN", "X_Longit", "Y_Latit"), 
           new = c("Mouse_ID", "Longitude", "Latitude"))
  
  HIJardaTable$Mouse_ID <- gsub(pattern = "SK", replacement = "SK_",x = HIJardaTable$Mouse_ID)
  
  # Let's remove the embryos (sadly, no interest for our parasitic studies...)
  HIJardaTable <- HIJardaTable[sapply(HIJardaTable$Mouse_ID, nchar) <= 7,]
  
  # Check if all rows are NA and delete these rows
  which(!rowSums(!is.na(HIJardaTable)))  
  
  # How many samples from Brandenburg do we have the HI for per year?
  table(HIJardaTable$Year)
  
  # CHECK DUPLICATE
  HIJardaTable$Mouse_ID[duplicated(HIJardaTable$Mouse_ID)]
  
  # Complement data with previous tables (2014, 2015) 
  diss2014.1 <- read.csv(paste0(pathToMyData, "HZ14_Mice 31-12-14_dissections.csv"),
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
  diss2014.2 <- read.csv(paste0(pathToMyData, "HZ14_Mice 31-12-14_genotypes.csv"),
                         stringsAsFactors=F)
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
  
  # Check if all rows are NA and delete these rows
  which(!rowSums(!is.na(mergedMiceTable)))  
  
  ## 2015 part 1
  diss2015.1 <- read.csv(paste0(pathToMyData, "Genotypes_Bav2015.csv"), 
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
  diss2015.2 <- read.csv(paste0(pathToMyData, "HZ15_Mice_Parasite.csv"), 
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
  
  # Check if all rows are NA and delete these rows
  which(!rowSums(!is.na(mergedMiceTable)))  
  
  ###################### 2016
  diss2016 <- read.csv(paste0(pathToMyData, "HZ16_Mice_18-07-16_dissections.csv"), 
                       na.strings=c(""," ","NA"), stringsAsFactors = F)[-c(1:2),]
  names(diss2016)[names(diss2016) %in% "ID_mouse"] <- "Mouse_ID"
  diss2016$Mouse_ID <- as.character(diss2016$Mouse_ID)
  
  ## Add worms 
  worms16 <- read.csv(paste0(pathToMyData, "HZ16_Worms.csv"), 
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
  
  # Check if all rows are NA and delete these rows
  which(!rowSums(!is.na(mergedMiceTable)))  
  
  mergedMiceTable <- fillGapsAfterMerge(mergedMiceTable)
  
  # Check if all rows are NA and delete these rows
  which(!rowSums(!is.na(mergedMiceTable)))  
  
  # add missing or wrong latitude/longitude
  loc2016 <- read.csv(paste0(pathToMyData, "Cleaned_HMHZ_2016_All.csv" ),
                      stringsAsFactors=F)
  loc2016 <- loc2016[names(loc2016) %in% c( "location", "GPS.coordinates.long", "GPS.coordinates.lat")]
  names(loc2016) <- c("Code", "longitude", "latitude")
  
  mergedMiceTable <- merge(mergedMiceTable, loc2016, by = "Code", all.x = TRUE)
  
  # 1. empty lon/lat
  mergedMiceTable$Longitude[is.na(mergedMiceTable$Longitude)] = 
    mergedMiceTable$longitude[is.na(mergedMiceTable$Longitude)]
  
  mergedMiceTable$Latitude[is.na(mergedMiceTable$Latitude)] = 
    mergedMiceTable$latitude[is.na(mergedMiceTable$Latitude)]
  
  # 2. wrong lon/lat
  mergedMiceTable$Latitude[!is.na(mergedMiceTable$Latitude) &
                             !is.na(mergedMiceTable$latitude) &
                             mergedMiceTable$Latitude != mergedMiceTable$latitude] =
    mergedMiceTable$latitude[!is.na(mergedMiceTable$Latitude) &
                               !is.na(mergedMiceTable$latitude) &
                               mergedMiceTable$Latitude != mergedMiceTable$latitude]
  
  
  mergedMiceTable$Longitude[!is.na(mergedMiceTable$Longitude) &
                              !is.na(mergedMiceTable$longitude) &
                              mergedMiceTable$Longitude != mergedMiceTable$longitude] =
    mergedMiceTable$longitude[!is.na(mergedMiceTable$Longitude) &
                                !is.na(mergedMiceTable$longitude) &
                                mergedMiceTable$Longitude != mergedMiceTable$longitude]
  # 3. delete duplicated rows
  mergedMiceTable <- unique(mergedMiceTable)
  
  ## **********************************************************
  ## 2017
  diss2017 <- read.csv(paste0(pathToMyData, "HZ17_September_Mice_Dissection.csv"),
                       na.strings=c(""," ","NA"), stringsAsFactors = F)
  
  # correction excel bullshit
  diss2017$Feces_weight <- as.numeric(as.character(diss2017$Feces_weight))
  
  diss2017$Feces_weight[diss2017$Feces_weight > 100 & !is.na(diss2017$Feces_weight)] <-
    diss2017$Feces_weight[diss2017$Feces_weight > 100 & !is.na(diss2017$Feces_weight)] / 1000
  
  names(diss2017)[names(diss2017) == "Ectoparasites"] <- "Flea"      
  
  ## Add worms 
  worms17 <- read.csv2(paste0(pathToMyData, "HZ17_September_Mice_Dissection_Jen_final.csv"),
                       stringsAsFactors = F)
  
  names(worms17)[names(worms17) %in% "Mesocestoides"] <- "Taenia_martis"
  
  ## merge worms and dissection table 2016
  diss2017 <- merge(diss2017, worms17, 
                    by = "Mouse_ID", all = TRUE)
  diss2017 <- fillGapsAfterMerge(diss2017)
  
  # merge
  mergedMiceTable <- merge(mergedMiceTable, diss2017, 
                           by = c("Mouse_ID"), all = T)
  
  # Check if all rows are NA and delete these rows
  which(!rowSums(!is.na(mergedMiceTable)))  
  
  mergedMiceTable <- fillGapsAfterMerge(mergedMiceTable)
  
  # Check if all rows are NA and delete these rows
  which(!rowSums(!is.na(mergedMiceTable)))  
  
  ## Uniformisation
  mergedMiceTable$Sex[grep("female*.", mergedMiceTable$Sex)] <- "F"
  mergedMiceTable$Sex[grep("male*.", mergedMiceTable$Sex)] <- "M"
  
  # Add old HI from previous jarda table
  oldJarda <- read.csv(paste0(pathToMyData, "HIforEH_May2017.csv"),
                       stringsAsFactors=F)
  
  # Uniformize IDs
  oldJarda$Mouse_ID <- oldJarda$PIN
  oldJarda$Mouse_ID <- gsub(pattern = "SK", replacement = "SK_",x = oldJarda$PIN)
  
  # Add ommited samples
  toadd <- oldJarda[!oldJarda$Mouse_ID %in% mergedMiceTable$Mouse_ID,]
  
  # Add samples without HI in mergedmicetable but with one in oldJarda
  missing <- mergedMiceTable$Mouse_ID[is.na(mergedMiceTable$HI)]
  
  toadd <- rbind(toadd, oldJarda[oldJarda$Mouse_ID %in% missing,])
  
  # mergedMiceTable[mergedMiceTable$Mouse_ID %in% missing 
  #                 & is.na(mergedMiceTable$Longitude) &
  #                   mergedMiceTable$Mouse_ID %in% toadd$Mouse_ID,]
  
  # setnames and merge
  setnames(toadd,
           old = c("BW", "L", "LCd", 
                   "Aspiculuris", "Syphacia", 
                   "Trichuris","Taenia"),
           new = c("Body_weight", "Body_length", "Tail_length", 
                   "Aspiculuris_tetraptera", "Syphacia_obvelata", 
                   "Trichuris_muris", "Taenia_taeniformis"))
  
  toadd = toadd[names(toadd)%in%c("REGion", "mtBamH", "Zfy2", "SRY1", "Y", "X332", "X347",
                                  "X65", "Tsx", "Btk", "Syap1", "Es1", "Gpd1", "Idh1",
                                  "Mpi", "Np", "Sod1", "Es1C", "Gpd1C", "Idh1C", "MpiC", "NpC",
                                  "Sod1C", "HI_NLoci", "HI", "Mouse_ID")]
  
  mergedMiceTable = merge(mergedMiceTable, toadd, by = "Mouse_ID", all = T)
  
  mergedMiceTable = fillGapsAfterMerge(mergedMiceTable)
  
  # correct error tail length
  mergedMiceTable$Tail_length <- as.numeric(mergedMiceTable$Tail_length)
  
  # delete duplicated rows
  mergedMiceTable <- unique(mergedMiceTable)
  
  duplicated = mergedMiceTable$Mouse_ID[duplicated(mergedMiceTable$Mouse_ID)]
  
  # Check if all rows are NA and delete these rows
  which(!rowSums(!is.na(mergedMiceTable)))  
  
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
  
  # Final cleaning, and save!
  mergedMiceTable$Longitude <- as.numeric(mergedMiceTable$Longitude)
  mergedMiceTable$Latitude <- as.numeric(mergedMiceTable$Latitude)
  
  # Check if all rows are NA and delete these rows
  which(!rowSums(!is.na(mergedMiceTable)))  
  
  ## Remove useless mice:
  
  # wildpark Schorfheide (not needed, test)
  wsh <- c(paste0("AA_000", 1:9), paste0("AA_00", 10:46))
  # apodemus caught in 2016
  apd <- c("A_0001", "A_0002", "A_0003")
  # useless info
  useless <- c(wsh, apd)
  
  mergedMiceTable <- mergedMiceTable[!(mergedMiceTable$Mouse_ID %in% useless),]
  
  # Check if all rows are NA and delete these rows
  which(!rowSums(!is.na(mergedMiceTable)))  
  
  # correct body length/weight
  mergedMiceTable$Body_length <- as.numeric(mergedMiceTable$Body_length)
  
  mergedMiceTable$Body_weight <- as.numeric(mergedMiceTable$Body_weight)
  
  # Check if all rows are NA and delete these rows
  which(!rowSums(!is.na(mergedMiceTable)))  
  
  # Manual correction
  mergedMiceTable$Body_weight[!is.na(mergedMiceTable$Body_weight) &
                                mergedMiceTable$Body_weight > 100] <-
    mergedMiceTable$Body_weight[!is.na(mergedMiceTable$Body_weight) &
                                  mergedMiceTable$Body_weight > 100] / 1000
  
  mergedMiceTable$Body_length[which(mergedMiceTable$Body_length < 20)] <- 
    mergedMiceTable$Body_length[which(mergedMiceTable$Body_length < 20)] * 10
  
  # Body condition index as log body mass/log body length (Hayes et al. 2014)
  mergedMiceTable$BCI <- log(mergedMiceTable$Body_weight) / log(mergedMiceTable$Body_length)
  
  # Correct wrong HI (>1)
  mergedMiceTable$HI[mergedMiceTable$HI > 1 & !is.na(mergedMiceTable$HI)] <- 
    mergedMiceTable$HI[mergedMiceTable$HI > 1 & !is.na(mergedMiceTable$HI)]/1000
  
  # Correct wrong Long/Lat: 1. wrong multiplicative factor
  mergedMiceTable$Longitude[mergedMiceTable$Longitude > 100 &
                              !is.na(mergedMiceTable$Longitude)] <- 
    mergedMiceTable$Longitude[mergedMiceTable$Longitude > 100 &
                                !is.na(mergedMiceTable$Longitude)] / 1000
  
  mergedMiceTable$Latitude[mergedMiceTable$Latitude > 100 & 
                             !is.na(mergedMiceTable$Latitude)] <- 
    mergedMiceTable$Latitude[mergedMiceTable$Latitude > 100 & 
                               !is.na(mergedMiceTable$Latitude)] / 1000
  
  # Correct wrong Long/Lat: 2. wrong lat/lon inversion
  mergedMiceTable$Longitude.temp <- mergedMiceTable$Longitude
  
  mergedMiceTable$Longitude[!is.na(mergedMiceTable$Longitude) &
                              mergedMiceTable$Longitude >= 30] =
    mergedMiceTable$Latitude[!is.na(mergedMiceTable$Longitude) &
                               mergedMiceTable$Longitude >= 30]
  
  mergedMiceTable$Latitude[!is.na(mergedMiceTable$Latitude) &
                             mergedMiceTable$Latitude <= 20] =
    mergedMiceTable$Longitude.temp[!is.na(mergedMiceTable$Latitude) &
                                     mergedMiceTable$Latitude <= 20]
  
  mergedMiceTable = 
    mergedMiceTable[ , -which(names(mergedMiceTable) %in% c("Longitude.temp"))]
  
  # add farm (TODO better localisation)
  mergedMiceTable$farm <- paste0(mergedMiceTable$Longitude, mergedMiceTable$Latitude)
  
  # Cluster by localities: rounded to about 700 meters ???...
  # mergedMiceTable$Latitude <- round(mergedMiceTable$Latitude, 2)
  # mergedMiceTable$Longitude <- round(mergedMiceTable$Longitude, 2)
  
  ## remove empty rows
  mergedMiceTable <- mergedMiceTable[!is.na(mergedMiceTable$Mouse_ID),]
  
  ## remove duplicated rows
  mergedMiceTable <- unique(mergedMiceTable)
  
  ## Correct duplicated mice by hand
  duplicated = mergedMiceTable$Mouse_ID[duplicated(mergedMiceTable$Mouse_ID)]
  
  ## 26 June 2018, add Jarda new csv
  missingMice = mergedMiceTable$Mouse_ID[is.na(mergedMiceTable$HI)]
  
  newCsv = read.csv(paste0(pathToMyData, "EmanuelData_26061018.csv"),
                    stringsAsFactors=F)
  newCsv$PIN = gsub("SK", "SK_", newCsv$PIN)
  
  ## Correct previous mistakes before merging for Latitude, Longitude, Transect, Sex, Code
  names(newCsv)[names(newCsv) %in% c("PIN", "X_Longit", "Y_Latit")] = 
    c("Mouse_ID", "Longitude", "Latitude")
  
  dataToAdd = newCsv[!newCsv$Mouse_ID %in% mergedMiceTable$Mouse_ID |
                       newCsv$Mouse_ID %in% missingMice,]
  
  mergedMiceTable =  merge(mergedMiceTable, dataToAdd, by = "Mouse_ID", all = T)
  
  mergedMiceTable = fillGapsAfterMerge(mergedMiceTable)
  
  # correct year manually
  mergedMiceTable$Year[
    mergedMiceTable$Mouse_ID %in% c("SK_2903", "SK_2904")] <- 2014
  mergedMiceTable$Year[
    mergedMiceTable$Mouse_ID %in% c("AA_0330", "AA_0450", "AA_0451", "AA_0452")] <- 2017
  
  # Add Eimeria information
  ## flotation
  flot <- read.csv("data/Eimeria_detection/FINALOocysts2015to2017.csv")
  
  ## how many neubauer cells were counted 
  flot$Ncells <- apply(flot[paste0("N_oocysts_sq", 1:8)], 1, function(x) sum(!is.na(x)))
  
  flot$OPG <- rowSums(flot[,paste0("N_oocysts_sq", 1:8)], na.rm = T) / flot$Ncells * 10000 /
    (flot$PBS_dil_in_mL * flot$Feces_g)
  
  # flot$Mouse_ID[!flot$Mouse_ID %in% miceTable$Mouse_ID]
  # SK_3174 only missing :(
  mergedMiceTable <- merge(mergedMiceTable, flot, all = T)
  
  ## qPCr
  qpcr <- read.csv("data/Eimeria_detection/FINALqpcrData_2016_2017_threshold5.csv")
  #qpcr$Mouse_ID[!qpcr$Mouse_ID %in% miceTable$Mouse_ID]
  #all there :)
  mergedMiceTable <- merge(mergedMiceTable, 
                           qpcr[c("Mouse_ID", "delta_ct_ilwe_MminusE", "delta_ct_cewe_MminusE", "observer")],
                           all = T)
  
  ## species identification
  species <- read.csv("data/Eimeria_detection/Eimeria_species_assignment_14_17.csv")
  names(species)[names(species) %in% "Species"] <- "eimeriaSpecies"
  
  # space error damn!
  species$Mouse_ID <- gsub(" ", "", as.character(species$Mouse_ID))
  mergedMiceTable <- merge(mergedMiceTable, 
                           species[c("Mouse_ID", "n18S_Seq", "COI_Seq", "ORF470_Seq", "eimeriaSpecies")],
                           by = "Mouse_ID", all = T)
  
  # clean Year
  mergedMiceTable[grep("A_00", mergedMiceTable$Mouse_ID),"Year"] <- 2016

  ## TO BE CHANGED WITH TISSUE RESULTS
  mergedMiceTable$eimeriaSpecies <- as.character(mergedMiceTable$eimeriaSpecies)
  mergedMiceTable$eimeriaSpecies[mergedMiceTable$Mouse_ID %in% "AA_0111"] <- "Double_ferrisi_vermiformis"
  mergedMiceTable$eimeriaSpecies[mergedMiceTable$Mouse_ID %in% "AA_0244"] <- "Double_tbd"
  mergedMiceTable$eimeriaSpecies[mergedMiceTable$Mouse_ID %in% "AA_0245"] <- "Double_tbd"
  mergedMiceTable$eimeriaSpecies[mergedMiceTable$Mouse_ID %in% "AA_0436"] <- "Double_tbd"
  mergedMiceTable$eimeriaSpecies[mergedMiceTable$Mouse_ID %in% "AA_0497"] <- "Double_tbd"
  mergedMiceTable$eimeriaSpecies <- as.factor(mergedMiceTable$eimeriaSpecies)
  
  return(mergedMiceTable)
}

# miceTable <- makeMiceTable()
# write.csv(miceTable, "data/MiceTable_fullEimeriaInfos_2014to2017.csv", row.names = F)
