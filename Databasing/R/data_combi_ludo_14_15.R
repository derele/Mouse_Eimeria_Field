raw.loc.2014 <- read.csv("HZ14_Mice 31-12-14_localities.csv")
raw.loc.2014 <- raw.loc.2014[, c(1:6, 13:15)]

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

raw.loc.2014$Latitude <- convert.deg(raw.loc.2014$Latitude)
raw.loc.2014$Longitude <- convert.deg(raw.loc.2014$Longitude)

genotypes.2014 <- read.csv("HZ14_Mice 31-12-14_genotypes.csv")[-1,]
X.col <- c("X332", "X347", "X65", "Tsx", "Btk",  "Syap1")


genotypes.2014$collapsed.GT <-
    apply(genotypes.2014[, X.col],
                         1, paste, collapse = "/")

get.HIX <- function (x){
    dom <- nchar(x) - nchar(gsub("d", "", x))
    mus <- nchar(x) - nchar(gsub("m", "", x))
    mus/(mus + dom)
}


genotypes.2014$HIX <-
    sapply(genotypes.2014$collapsed.GT, get.HIX)

loc.HIX.2014 <- get.HIX(
    tapply(genotypes.2014$collapsed.GT, genotypes.2014$Code,
           function (x){
               paste(x, collapse="/")}))

loc.2014 <- merge(raw.loc.2014, loc.HIX.2014, by.x="Code", by.y = 0)


ludo <- read.csv("dureje_et_al._supplementary_table_s1.csv",
                 skip=1)

ludo$Latitude <- convert.deg(ludo$Latitude)
ludo$Longitude <- convert.deg(ludo$Longitude)

ludo <- ludo[!ludo$Longitude==0, ]
ludo <- ludo[!ludo$Latitude==0, ]

library(ggmap)

area <- get_map(location =
                    c(min(ludo$Longitude),
                      min(ludo$Latitude),
                      max(ludo$Longitude),
                      max(ludo$Latitude)),
                source = "stamen", maptype="toner-lite")

png("figures/Ludo.png", units = "in", width = 6, height = 7, res=300)
ggmap(area,  zoom = 15) +
    geom_point(data = subset(ludo, !is.na(ludo$HIX)),
               aes(Longitude, Latitude, color=HIX), size=2, alpha=0.5) +
                   scale_color_gradient("Hybrid\nindex", high="red",low="blue")
dev.off()


png("figures/Ludo_plus_Berlin.png", units = "in", width = 6,
    height = 7, res=300)
ggmap(area,  zoom = 15) +
    geom_point(data = subset(ludo, !is.na(ludo$HIX)),
               aes(Longitude, Latitude, color=HIX),size=2, alpha=0.5) +
        scale_color_gradient("Hybrid\nindex", high="red",low="blue") +
            geom_point(data = subset(loc.2014, !is.na(loc.2014$y)),
                       aes(Longitude, Latitude, color=y), size=2, alpha=0.5) 
dev.off()


area.B <- get_map(location =
                      c(12.4,
                        52,
                        max(loc.2014$Longitude),
                        max(loc.2014$Latitude)),
                source = "stamen", maptype="toner-lite")

png("figures/Berlin.png", units="in", width = 6, height = 6,
    res=300)
ggmap(area.B,  zoom = 15) +
            geom_point(data = subset(loc.2014, !is.na(loc.2014$y)),
                       aes(Longitude, Latitude, color=y), size=2) +
                scale_color_gradient("Hybrid\nindex", high="red",low="blue") 
dev.off()

raw.loc.2015 <- read.csv("HZ15_Mice_localities.csv")

raw.loc.2015$X_Map <- as.numeric(as.character(raw.loc.2015$X_Map))
raw.loc.2015$Y_Map <- as.numeric(as.character(raw.loc.2015$Y_Map))



png("figures/Ludo_plus_Berlin_plus_Bavaria.png", units="in", width = 6, height = 7,  res=300)
ggmap(area,  zoom = 15) +
    geom_point(data = raw.loc.2015,
               aes(X_Map, Y_Map), color="green", size=2)+
        geom_point(data = subset(ludo, !is.na(ludo$HIX)),
                   aes(Longitude, Latitude, color=HIX),size=2, alpha=0.5) +
                geom_point(data = subset(loc.2014, !is.na(loc.2014$y)),
                           aes(Longitude, Latitude, color=y), size=2, alpha=0.5)+
                    scale_color_gradient("Hybrid\nindex", high="red",low="blue")
dev.off()


HZ15G <- read.csv("HZ15_Mice_genotypes.csv")

HZ15G$collapsed.GT <-
    apply(HZ15G[, X.col],
          1, paste, collapse = "/")


HZ15G$HIX <- sapply(HZ15G$collapsed.GT, get.HIX)


raw.loc.2015 <- read.csv("HZ15_Mice_localities.csv")

raw.loc.2015$X_Map <- as.numeric(as.character(raw.loc.2015$X_Map))
raw.loc.2015$Y_Map <- as.numeric(as.character(raw.loc.2015$Y_Map))


loc.HIX.2015 <- get.HIX(
    tapply(HZ15G$collapsed.GT, HZ15G$Code,
           function (x){
               paste(x, collapse="/")}))

loc.2015 <- merge(raw.loc.2015, loc.HIX.2015, by.x="Code", by.y = 0)
names(loc.2015)[names(loc.2015)%in%"y"] <- "HIX"

HZ15G.Bav <- read.csv("Genotypes_Bav2015.csv")

HZ15G.Bav$collapsed.GT <-
    apply(HZ15G.Bav[, X.col],
          1, paste, collapse = "/")

HZ15G.Bav$HIX <- sapply(HZ15G.Bav$collapsed.GT, get.HIX)

loc.HIX.2015.Bav <- get.HIX(
    tapply(HZ15G.Bav$collapsed.GT, HZ15G.Bav$Code,
           function (x){
               paste(x, collapse="/")}))

loc.2015 <- merge(raw.loc.2015, loc.HIX.2015.Bav, by.x="Code", by.y = 0)
names(loc.2015)[names(loc.2015)%in%"y"] <- "HIX"

png("figures/Ludo_plus_Berlin15_plus_Bavaria.png", units="in", width = 6, height = 7,  res=300)
ggmap(area,  zoom = 15) +
        geom_point(data = subset(ludo, !is.na(ludo$HIX)),
                   aes(Longitude, Latitude, color=HIX),size=2, alpha=0.5) +
                geom_point(data = subset(loc.2014, !is.na(loc.2014$y)),
                           aes(Longitude, Latitude, color=y), size=2, alpha=0.5)+
                    scale_color_gradient("Hybrid\nindex", high="red",low="blue") +
                        geom_point(data = loc.2015,
                                   aes(X_Map, Y_Map, color=HIX), size=2, alpha=0.5) 
dev.off()


png("figures/Berlin_2015.png", units="in", width = 6, height = 6,  res=300)
ggmap(area.B,  zoom = 15) +
    scale_color_gradient("Hybrid\nindex", high="red",low="blue")+
        geom_point(data = loc.2015,
                   aes(X_Map, Y_Map, color=HIX, alpha=0.5))+
            geom_point(data = subset(loc.2014, !is.na(loc.2014$y)),
                       aes(Longitude, Latitude, color=y), alpha=0.5)
dev.off()

png("figures/Bavaria_2015.png", units="in", width = 6, height = 6,  res=300)
ggmap(area.B,  zoom = 15) +
    scale_color_gradient("Hybrid\nindex", high="red",low="blue")+
        geom_point(data = loc.2015,
                   aes(X_Map, Y_Map, color=HIX, alpha=0.5))+
            geom_point(data = subset(loc.2014, !is.na(loc.2014$y)),
                       aes(Longitude, Latitude, color=y), alpha=0.5)
dev.off()

png("figures/Berlin_2015_labeled.png", units="in", width = 6, height = 6,  res=300)
ggmap(area.B,  zoom = 15) +
    scale_color_gradient("Hybrid\nindex", high="red",low="blue")+
        geom_point(data = loc.2015,
                   aes(X_Map, Y_Map, color=HIX, alpha=0.5))+
            geom_point(data = subset(loc.2014, !is.na(loc.2014$y)),
                       aes(Longitude, Latitude, color=y), alpha=0.5) +
                geom_text(data = loc.2015,
                           aes(X_Map, Y_Map, label=Code), size = 1)
dev.off()





