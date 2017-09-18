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

## rodents caught
sum(Trap$Number_traps_set)

## localities set
nrow(Trap[!duplicated(Trap[, c("Latitude", "Longitude")]),
          c("Latitude", "Longitude")])

## localities caught
nrow(Trap[Trap$Number_rodents_caught>0 &
              !duplicated(Trap[, c("Longitude", "Latitude")]),])

nrow(Trap[Trap$Number_mus_caught>0 &
              Trap$Number_rodents_caught > Trap$Number_mus_caught &
                  !duplicated(Trap[, c("Longitude", "Latitude")]),])

## distribution of trapping success
pdf("/home/ele/Dist_mus_caught.pdf")
ggplot(Trap, aes(Number_mus_caught)) +
    geom_histogram(binwidth=1)
dev.off()

## distribution of trapping attempts
pdf("/home/ele/Dist_traps_set.pdf")
ggplot(subset(Trap, !Team%in%"X"),
       aes(Number_traps_set, ..density.., color=Team)) +
#    geom_histogram(binwidth=1) +
        geom_density()
dev.off()

pdf("/home/ele/Traps_vs_Rodents.pdf")
ggplot(Trap, aes(Number_traps_set, Number_rodents_caught,
                 color=Team)) +
                     geom_jitter(width=0.2) +
                         stat_smooth(method="lm", se=FALSE)
dev.off()

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

People.tab$efficiency <- People.tab$nMus/People.tab$nTraps

People.tab[order(People.tab$nTraps), ]
People.tab[order(People.tab$nMus), ]
People.tab[order(People.tab$nRodents), ]
People.tab[order(People.tab$nLoc), ]
People.tab[order(People.tab$nLocM), ]
People.tab[order(People.tab$efficiency), ]

Farm_eco <- read.csv("/home/ele/Farm_ecology.csv")
Farm_eco$Longitude <- gsub( ",", ".", Farm_eco$Longitude)
Farm_eco$Latitude <- gsub( ",", ".", Farm_eco$Latitude)

Farm_table <- merge(Trap, Farm_eco)

write.csv(Farm_table, "/home/ele/Merged_Farm.csv")

apply(Farm_table[, 17:ncol(Farm_table)], 2, table)

pdf("/home/ele/Getreidefeld.pdf")
ggplot(Farm_table, aes(as.factor(Feld.Getreide.), Number_mus_caught)) +
    geom_boxplot() + 
    geom_jitter(width=0.2)
dev.off()

wilcox.test(Farm_table$Number_mus_caught[Farm_table$Feld.Getreide.==0],
            Farm_table$Number_mus_caught[Farm_table$Feld.Getreide.>0])

pdf("/home/ele/Bebautes_Land.pdf")
ggplot(Farm_table, aes(as.factor(Bebautes.Land), Number_mus_caught)) +
    geom_boxplot() + 
    geom_jitter(width=0.2)
dev.off()

wilcox.test(Farm_table$Number_mus_caught[Farm_table$Bebautes.Land<5],
            Farm_table$Number_mus_caught[Farm_table$Bebautes.Land>=5])

pdf("/home/ele/Wald.pdf")
ggplot(Farm_table, aes(as.factor(Wald), Number_mus_caught)) +
    geom_boxplot() + 
    geom_jitter(width=0.2)
dev.off()

wilcox.test(Farm_table$Number_mus_caught[Farm_table$Wald<5],
            Farm_table$Number_mus_caught[Farm_table$Wald>=5])


ggplot(Farm_table, aes(as.factor(Wald), Number_rodents_caught-Number_mus_caught)) +
    geom_boxplot() + 
    geom_jitter(width=0.2)

wilcox.test(Farm_table$Number_mus_caught[Farm_table$Wald<5],
            Farm_table$Number_mus_caught[Farm_table$Wald>=5])

pdf("/home/ele/Ratten.pdf")
ggplot(Farm_table, aes(as.factor(Ratten), Number_mus_caught)) + geom_boxplot() + geom_jitter()
dev.off()

wilcox.test(Farm_table$Number_mus_caught[Farm_table$Ratten<1],
            Farm_table$Number_mus_caught[Farm_table$Ratten>0])

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
