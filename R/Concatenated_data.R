## Concatenate data 
setwd("/home/alice/git_projects/Mouse_Eimeria_Databasing/raw_data/")

######### Functions to be used:
## Convert GPS coordinate in decimals format
convert.deg <-function(c){
    z <- lapply(strsplit(as.character(c),
                         "° *|' *|(\" *|$)"), as.numeric)
    ## fill to seconds if absent
    zz <- lapply(z, function(X) {
                     c(X, rep(0, times = 3 - length(X)))
                 })
    dec <- lapply(zz, function (x) x[1] + x[2]/60 + x[3]/3600)
    return(unlist(dec))
} 
## Calculate the hybrid index
get.HIX <- function (x){
    dom <- nchar(x) - nchar(gsub("d", "", x))
    mus <- nchar(x) - nchar(gsub("m", "", x))
    mus/(mus + dom)
}
## map rows from different data.frames by their common column
rbind.match.columns <- function(input1, input2) {
    n.input1 <- ncol(input1)
    n.input2 <- ncol(input2)
    if (n.input2 < n.input1) {
        TF.names <- which(names(input2) %in% names(input1))
        column.names <- names(input2[, TF.names])
    } else {
        TF.names <- which(names(input1) %in% names(input2))
        column.names <- names(input1[, TF.names])
    }
    return(rbind(input1[, column.names], input2[, column.names]))
}
#########

#################
## Data from 2014
## locations
raw.loc.2014 <- read.csv("HZ14_Mice 31-12-14_localities.csv")
raw.loc.2014 <- raw.loc.2014[, c(1:6, 13:15)]
raw.loc.2014$Latitude <- convert.deg(raw.loc.2014$Latitude)
raw.loc.2014$Longitude <- convert.deg(raw.loc.2014$Longitude)

## genotypes
genotypes.2014 <- read.csv("HZ14_Mice 31-12-14_genotypes.csv")[-1,]
X.col <- c("X332", "X347", "X65", "Tsx", "Btk",  "Syap1")
## collapse m or d alleles in 1 column
genotypes.2014$collapsed.GT <-
    apply(genotypes.2014[, X.col],
                         1, paste, collapse = "/")
## calculate hybrid index per mouse
genotypes.2014$HIX <-
    sapply(genotypes.2014$collapsed.GT, get.HIX)
## calculate hybrid index per locality
loc.HIX.2014 <- get.HIX(
    tapply(genotypes.2014$collapsed.GT, genotypes.2014$Code,
           function (x){
               paste(x, collapse="/")}))

loc.2014 <- merge(raw.loc.2014, loc.HIX.2014, by.x="Code", by.y = 0)

## dissection data
dissection.2014 <- read.csv("HZ14_Mice 31-12-14_dissections.csv")

#################
## Data from Ludo
ludo <- read.csv("dureje_et_al._supplementary_table_s1.csv",
                 skip=1)

ludo$Latitude <- convert.deg(ludo$Latitude)
ludo$Longitude <- convert.deg(ludo$Longitude)

ludo <- ludo[!ludo$Longitude==0, ]
ludo <- ludo[!ludo$Latitude==0, ]

#################
## Data from 2015
## locations
raw.loc.2015 <- read.csv("HZ15_Mice_localities.csv")
raw.loc.2015$X_Map <- as.numeric(as.character(raw.loc.2015$X_Map))
raw.loc.2015$Y_Map <- as.numeric(as.character(raw.loc.2015$Y_Map))

## genotypes part 1
HZ15G <- read.csv("HZ15_Mice_genotypes.csv")
## collapse m or d alleles in 1 column
HZ15G$collapsed.GT <-
    apply(HZ15G[, X.col], ## apply to the 6 markers column
          1, paste, collapse = "/")
## calculate hybrid index per mouse
HZ15G$HIX <- sapply(HZ15G$collapsed.GT, get.HIX)
## calculate hybrid index per locality
loc.HIX.2015 <- get.HIX(
    tapply(HZ15G$collapsed.GT, HZ15G$Code,
           function (x){
               paste(x, collapse="/")}))
loc.2015 <- merge(raw.loc.2015, loc.HIX.2015, by.x="Code", by.y = 0)
names(loc.2015)[names(loc.2015)%in%"y"] <- "HIX"

## genotypes part 2 (Bavaria)
HZ15G.Bav <- read.csv("Genotypes_Bav2015.csv")
## collapse m or d alleles in 1 column
HZ15G.Bav$collapsed.GT <-
    apply(HZ15G.Bav[, X.col],
          1, paste, collapse = "/")
## calculate hybrid index per mouse
HZ15G.Bav$HIX <- sapply(HZ15G.Bav$collapsed.GT, get.HIX)
## calculate hybrid index per locality
loc.HIX.2015.Bav <- get.HIX(
    tapply(HZ15G.Bav$collapsed.GT, HZ15G.Bav$Code,
           function (x){
               paste(x, collapse="/")}))
loc.2015 <- merge(raw.loc.2015, loc.HIX.2015.Bav, by.x="Code", by.y = 0)
names(loc.2015)[names(loc.2015)%in%"y"] <- "HIX"

## dissection data
dissection.2015 <- read.csv("HZ14_Mice 31-12-14_dissections.csv")
parasites.2015 <- read.csv("HZ15_Mice_Parasite.csv")

#################
## Data from 2016
## locations
## table total (even no mouse captured)
extended.2016 <- read.csv("Cleaned_HMHZ_2016_All.csv")
## table with mice captured
raw.2016 <-read.csv("Cleaned_HMHZ_2016_MmOnly.csv") 
raw.2016$HZ <- "HZ_BR"
raw.2016 <- raw.2016[-c(1,2,7,8,9,10,11,12,13,14,15)]
raw.loc.2016 <- unique(raw.2016)

## genotypes
## .... to be completed

## dissection data
dissection.2016 <- read.csv("HZ16_Mice_18-07-16_dissections.csv")

###################
## Concatenated Map
library(ggmap)

map_2014 <- loc.2014
map_2014$Year <- "2014"
map_2014$HIX <- NA

map_ludo <- ludo
map_ludo$Year <- "compiled data"
map_ludo$Code <- NA

map_2015 <- loc.2015
names <- names(map_2015)
names(map_2015) <- c(names[1:6],"Longitude", "Latitude", names[9:25])
map_2015$Year <- "2015"

map_2016 <- raw.loc.2016
names(map_2016) <- c("Code", "Address", "Latitude", "Longitude", "HZ")
map_2016$Year <- "2016"
map_2016$HIX <- NA

a <- rbind.match.columns(map_2014, map_2015)
b <- rbind.match.columns(a, map_ludo)
map_all <- rbind.match.columns(b, map_2016)

## one general map
area <- get_map(location =
                    c(min(map_all$Longitude),
                      min(map_all$Latitude),
                      max(map_all$Longitude),
                      max(map_all$Latitude)),
                source = "stamen", maptype="toner-lite")

## install.packages("ggmap", type = "source")
png("Concatenated_map_feb2017.png", units = "in", width = 6,
    height = 7, res=300)
ggmap(area,  zoom = 16)+
    geom_point(data = map_all,
               aes(Longitude, Latitude, color=HIX, pch=Year), size=5, alpha=0.7) +
    scale_color_gradient("Hybrid\nindex", high="red",low="blue")
dev.off()                         
