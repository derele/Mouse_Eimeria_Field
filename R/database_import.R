## Functions we need to run this
source("R/database_import_functions.R")


######## Locations ###############
raw.loc.2014 <- read.csv("HZ14_Mice 31-12-14_localities.csv")
raw.loc.2014$Latitude <- convert.deg(raw.loc.2014$Latitude)
raw.loc.2014$Longitude <- convert.deg(raw.loc.2014$Longitude)

raw.loc.2015 <- read.csv("HZ15_Mice_localities.csv")
raw.loc.2015$Longitude <- as.numeric(as.character(raw.loc.2015$X_Map))
raw.loc.2015$Latitude <- as.numeric(as.character(raw.loc.2015$Y_Map))


store <- c("Code", "Continent", "State", "Country",
           "Locality", "Latitude", "Longitude")

loc.raw <- rbind(raw.loc.2014[, store], raw.loc.2015[, store])

loc.consistent <- do.call(rbind, by(loc.raw, loc.raw$Code, is.geo.consistent))

loc.inconsistent <- do.call(rbind, by(loc.raw, loc.raw$Code,
                                      is.geo.consistent, consistent=FALSE))

loc.consistent <- rbind(loc.consistent, loc.inconsistent[2,])

loc <- loc.consistent[!duplicated(loc.consistent$Code),]

########## DB setup #######################

library(DBI)
library(RSQLite)

## Push to SQLite
db <- dbConnect(SQLite(), dbname = "/SAN/db/Test.sqlite")
dbWriteTable(conn = db, name = "Locations", value = loc, row.names = FALSE)

dbListFields(db, "Locations")

foo <- dbReadTable(db, "Locations")

########## Mouse #########################
GT.2014 <- read.csv("HZ14_Mice 31-12-14_genotypes.csv")
DS.2014 <- read.csv("HZ14_Mice 31-12-14_dissections.csv")

## Needed_not <- read.csv("HZ15_Mice_dissections.csv")
GT.2015 <- read.csv("Genotypes_Bav2015.csv")
DS.2015 <- read.csv("HZ15_Mice_Parasite.csv")


