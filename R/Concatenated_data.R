## Concatenate data Feb 2017
rm(list=ls())

setwd("/home/alice/git_projects/Mouse_Eimeria_Databasing/")

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

####################
## Cleaned locations
cleaned_loc <- read.csv("output_data/all_clean_localities.csv")
cleaned_loc <- cleaned_loc[-1]

## Add dereje data
ludo <- read.csv("raw_data/dureje_et_al._supplementary_table_s1.csv",
                 skip=1)
ludo$Latitude <- convert.deg(ludo$Latitude)
ludo$Longitude <- convert.deg(ludo$Longitude)
ludo <- ludo[!ludo$Longitude==0, ]
ludo <- ludo[!ludo$Latitude==0, ]
ludo$Loc. <- paste("Loc",ludo$Loc., sep="_")
names(ludo)[names(ludo)%in%"Loc."] <- "Code"

All_loc <- rbind.match.columns(cleaned_loc, ludo)

#############
### Genotypes
genotypes.2014 <- read.csv("raw_data/HZ14_Mice 31-12-14_genotypes.csv")[-1,]
X.col <- c("X332", "X347", "X65", "Tsx", "Btk",  "Syap1")
genotypes.2014$PIN <- paste0(genotypes.2014$ID,genotypes.2014$PIN)

genotypes.2015.1 <- read.csv("raw_data/HZ15_Mice_genotypes.csv")
genotypes.2015.1$PIN <- paste0(genotypes.2015.1$ID,genotypes.2015.1$PIN)
genotypes.2015.1$Year <- "2015"

genotypes.2015.2 <- read.csv("raw_data/Genotypes_Bav2015.csv")

## genotypes.2016 <- to be continued...

MyGen <- function(x) {
    cbind(x[,X.col], x$PIN,x$Code, x$Year)
}

Gen <- rbind(MyGen(genotypes.2014),
             MyGen(genotypes.2015.1),
             MyGen(genotypes.2015.2))

## collapse m or d alleles in 1 column
Gen$collapsed.GT <-
    apply(Gen[, X.col],1, paste, collapse = "/")

## calculate hybrid index per mouse
Gen$HIX <-
    sapply(Gen$collapsed.GT, get.HIX)

names(Gen)[c(7,8,9)] <- c("PIN", "Code", "Year")

## calculate hybrid index per locality
loc.HIX <- get.HIX(
    tapply(Gen$collapsed.GT, Gen$Code,
           function (x){
               paste(x, collapse="/")}))

## Add dereje data
ludo$Year <- "review derele"

All_Gen <- rbind.match.columns(Gen, ludo)

#######
## To remove when genotype 2016 will arrive
## to use before getting the genotypes data :
dissection.2016 <- read.csv("raw_data/HZ16_Mice_18-07-16_dissections.csv")
prev.2016 <- data.frame(unique(dissection.2016$Code), "2016", NA)
names(prev.2016) <- names(All_Gen)
All_Gen <- rbind(All_Gen, prev.2016)
#######

## remove duplicates
All_Gen <- unique(All_Gen)

Total_data <- merge(All_Gen, All_loc)

###################
## Concatenated Map
library(ggmap)

## one general map
area <- get_map(location =
                    c(min(Total_data$Longitude),
                      min(Total_data$Latitude),
                      max(Total_data$Longitude),
                      max(Total_data$Latitude)),
                source = "stamen", maptype="toner-lite")

## install.packages("ggmap", type = "source")
#png("/home/alice/git_projects/Mouse_Eimeria_Databasing/output_data/Concatenated_map_feb2017.png", units = "in", width = 6, height = 7, res=300)
ggmap(area,  zoom = 16)+
    geom_point(data = Total_data,
               aes(Longitude, Latitude, color=HIX, pch=Year), size=3, alpha=0.5) +
    scale_color_gradient("Hybrid\nindex", high="red",low="blue")
#dev.off()                         
