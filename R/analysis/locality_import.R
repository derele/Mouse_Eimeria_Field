## Functions we need to run this
source("R/database_import_functions.R")


######## Locations ###############
raw.loc.2014 <- read.csv("raw_data/HZ14_Mice 31-12-14_localities.csv")
raw.loc.2014$Latitude <- convert.deg(raw.loc.2014$Latitude)
raw.loc.2014$Longitude <- convert.deg(raw.loc.2014$Longitude)

raw.loc.2015 <- read.csv("raw_data/HZ15_Mice_localities.csv")
raw.loc.2015$Longitude <- as.numeric(as.character(raw.loc.2015$X_Map))
raw.loc.2015$Latitude <- as.numeric(as.character(raw.loc.2015$Y_Map))

raw.loc.2016 <- read.csv("raw_data/Cleaned_HMHZ_2016_MmOnly.csv")
raw.loc.2016$Latitude <- as.numeric(as.character(raw.loc.2016$GPS.coordinates.long))
raw.loc.2016$Longitude <- as.numeric(as.character(raw.loc.2016$GPS.coordinates.lat))

names(raw.loc.2016)[names(raw.loc.2016)%in%"location"] <- "Code"

store <- c("Code", "Latitude", "Longitude")

loc.raw <- rbind(raw.loc.2014[, store],
                 raw.loc.2015[, store],
                 raw.loc.2016[, store])

### Correct a bavarian locality which was labeled like one in Brandenburg
loc.raw[loc.raw$Code%in%"STEG1" &
        loc.raw$Latitude%in%"49.8825" &
        loc.raw$Longitude%in%"11.72"
       , "Code"] <- "SEYB"

### relable a minor locality with different buildings of the sampe property:
loc.raw[loc.raw$Code%in%"SCHEN8h" &
        loc.raw$Latitude%in%"52.273117" &
        loc.raw$Longitude%in%"13.599425" 
      , "Code"] <- "SCHEN8f"

### Localities sampled in 2015 by Ludo and Oli in Poland:


levels(loc.raw$Code) <- c(levels(loc.raw$Code), paste0("LudPol", 1:3))
loc.raw[loc.raw$Code%in%"", "Code"] <- paste0("LudPol", 1:3)

loc.consistent <- do.call(rbind, by(loc.raw, loc.raw$Code, is.geo.consistent,
                                    granularity=0.001))

loc.inconsistent <- do.call(rbind, by(loc.raw, loc.raw$Code,
                                      is.geo.consistent, granularity=0.001,
                                      consistent=FALSE))

loc <- loc.consistent[!duplicated(loc.consistent$Code),]

write.csv(loc, "output_data/all_clean_localities.csv")





########## DB setup #######################

## library(DBI)
## library(RSQLite)

## ## Push to SQLite
## db <- dbConnect(SQLite(), dbname = "/SAN/db/Test.sqlite")
## dbWriteTable(conn = db, name = "Locations", value = loc, row.names = FALSE)

## dbListFields(db, "Locations")

## foo <- dbReadTable(db, "Locations")

## ########## Mouse #########################
## GT.2014 <- read.csv("HZ14_Mice 31-12-14_genotypes.csv")
## DS.2014 <- read.csv("HZ14_Mice 31-12-14_dissections.csv")

## ## Needed_not <- read.csv("HZ15_Mice_dissections.csv")
## GT.2015 <- read.csv("Genotypes_Bav2015.csv")
## DS.2015 <- read.csv("HZ15_Mice_Parasite.csv")


