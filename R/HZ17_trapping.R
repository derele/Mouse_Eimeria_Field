library(plyr)
library(ggplot2)

## A trapping table with our new data
Trap <- read.csv("raw_data/HZ17_Mice_Trap.csv")
table(Trap$Number_mus_caught > 0)
Trap$Is_mus_caught <- Trap$Number_mus_caught > 0

## check double sampled localities 
dup.lat <- Trap[duplicated(Trap[, c("Latitude")]), c("Latitude")]
dup.long <- Trap[duplicated(Trap[, c("Longitude")]), c("Longitude")]  
Trap[Trap$Longitude%in%dup.long, ]

## check Mus are a subset of rodents
table(Trap$Number_rodents_caught >= Trap$Number_mus_caught)

## A locality table (trapping only tracking positive data) for pre 2017
Mice_till_17 <- read.csv("output_data/genDF_august2017.csv")

Loc_till_17 <- ddply (Mice_till_17, .(Longitude, Latitude),
   function(x){
       cbind(Year=unique(x$Year), HI=mean(x$HI), n.mice= nrow(x))
   })

## For use on interactive map
write.csv(Trap, "raw_data/HZ17_Mice_Trap_4map.csv")
write.csv(Loc_till_17, "output_data/HZ14-16_localities.csv")

## house mice caught
sum(Trap$Number_mus_caught)

## rodents caught
sum(Trap$Number_rodents_caught)

## localities set
nrow(Trap[!duplicated(Trap[, c("Latitude", "Longitude")]),
          c("Latitude", "Longitude")])

## localities caught
nrow(Trap[Trap$Number_mus_caught>0 &
              !duplicated(Trap[, c("Longitude", "Latitude")]),])

nrow(Trap[Trap$Number_mus_caught>0 &
              Trap$Number_rodents_caught > Trap$Number_mus_caught &
                  !duplicated(Trap[, c("Longitude", "Latitude")]),])

People <- strsplit(as.character(Trap$People), ", ?", perl=TRUE)

u.People <- unique(unlist(People))

People.tab <- t(sapply(u.People, function (name){
    select <- unlist(lapply(People, function (x) name %in% x))
    Tselect <- Trap[select, ]
    Ntraps <- sum(Tselect$Number_traps_set)
    Nrod <- sum(Tselect$Number_rodents_caught)
    Nmus <- sum(Tselect$Number_mus_caught)
    Nloc <- nrow(Tselect[!duplicated(Tselect[, c("Longitude", "Latitude")]),])
    NlocM <- nrow(Tselect[Tselect$Number_mus_caught>0 &
                          !duplicated(Tselect[, c("Longitude",
                                                  "Latitude")]),])
    cbind(Ntraps, Nrod, Nmus, Nloc, NlocM)
}))

People.tab <- as.data.frame(People.tab)

colnames(People.tab) <- c("nTraps", "nRodents", "nMus", "nLoc", "nLocM")

People.tab[order(People.tab$nTraps), ]
People.tab[order(People.tab$nMus), ]
People.tab[order(People.tab$nRodents), ]
People.tab[order(People.tab$nLoc), ]
People.tab[order(People.tab$nLocM), ]

is.loc.in.table <- function(Table,
                            coordinates,
                            granularity=0.001){
    apply(Table, 1, function(x){
        diff <- as.numeric(as.character(x[c("Longitude", "Latitude")])) -
            coordinates
        all(abs(diff) < granularity)})
}


all.net <- apply(Loc_till_17, 1, function (x) {
     all.net <- is.loc.in.table(Trap,
                                x[c("Longitude", "Latitude")])
     all.net
})

r.new <- apply(r.new, 1, paste, collapse="|")

r.old <- round(unique(Loc_till_17[, c("Longitude", "Latitude")]), 3)
r.old <- apply(r.old, 1, paste, collapse="|")


table(r.new%in%r.old)
