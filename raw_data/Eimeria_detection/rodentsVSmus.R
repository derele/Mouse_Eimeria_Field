rodentTable <- read.csv("../../Data_important/MiceTable_2014to2017.csv")
PCRmice <- read.csv("../../Mouse_Eimeria_Databasing/raw_data/Eimeria_detection/Inventory_contents_all.csv")
PCRotherRodents <- read.csv("../../Mouse_Eimeria_Databasing/raw_data/Eimeria_detection/Other_rodents/rawdata_other_rodents.csv", 
                            +                             na.strings = "")

# Which other rodents are infected?
idOR <- PCRotherRodents$Mouse_ID[!is.na(PCRotherRodents$Eimeria_ident)]
lat <- PCRotherRodents$Latitude[!is.na(PCRotherRodents$Eimeria_ident)]
lon <- PCRotherRodents$Longitude[!is.na(PCRotherRodents$Eimeria_ident)]
host <- PCRotherRodents$Host[!is.na(PCRotherRodents$Eimeria_ident)]

ORdf <- data.frame(idOR, lat, lon, host)

# Which mus are infected?
df <- PCRmice[PCRmice$X13_18S_Seq %in% "positive" |
                +                    PCRmice$X14_COI_Seq %in% "positive" |
                +                    PCRmice$X16_ORF470_Seq %in% "positive",]
idMUS <- df$X3_ID_mouse
idMUS <- gsub(" ", "", idMUS)
idMUS <- unique(idMUS)

# Lat/lon for these mice?
MUSdf <- rodentTable[rodentTable$Mouse_ID %in% idMUS, c("Mouse_ID", "Latitude", "Longitude")]
MUSdf <- MUSdf[!is.na(MUSdf$Latitude) & !is.na(MUSdf$Longitude), ]

# formatting before merge
MUSdf$lon <- MUSdf$Longitude
MUSdf$lat <- MUSdf$Latitude

# MERGE AND SEEEEE
OMG <- merge(MUSdf, ORdf, all = T)

OMG$species[is.na(OMG$Latitude)] <- as.character(OMG$host[is.na(OMG$Latitude)])
OMG$species[!is.na(OMG$Latitude)] <- "mus"

library(geosphere)
musmat <- as.matrix(data.frame(MUSdf$lat, MUSdf$lon))
ormat <- as.matrix(data.frame(ORdf$lat, ORdf$lon))

distanceMatrix <- distm(musmat, ormat, fun = distHaversine)/1000 # in km
distances <- data.frame(distanceMatrix)
names(distances) <- ORdf$idOR
distances$musid <- MUSdf$Mouse_ID


# all distance rodent
allmustable <- rodentTable[grep("AA", rodentTable$Mouse_ID),]

allmusmat <- as.matrix(data.frame(allmustable$Latitude, allmustable$Longitude))
distanceMatrix2 <- distm(allmusmat, ormat, fun = distHaversine)/1000 # in km

distances2 <- data.frame(distanceMatrix2)
names(distances2) <- ORdf$idOR
distances2$musid <- allmustable$Mouse_ID

distances2 <- distances2[!is.na(distances2$ZZ_0021),]

distances2[distances2[,1] < 1,]
distances2[distances2[,2] < 1,]
distances2[distances2[,3] < 1,]
distances2[distances2[,4] < 1,]
distances2[distances2[,5] < 1,]
distances2[distances2[,6] < 1,]
distances2[distances2[,7] < 1,]
distances2[distances2[,8] < 1,]
