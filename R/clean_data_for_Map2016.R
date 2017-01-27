library(ggplot2)
library(ggmap)

messy <- read.csv("HMHZ_2016_tofillaftercatching.csv", sep=",")

clean <- messy

clean[, "GPS.coordinates.long"] <- gsub(",", ".", clean[, "GPS.coordinates.long"])
clean[, "GPS.coordinates.lat"] <- gsub(",", ".", clean[, "GPS.coordinates.lat"])

clean <- as.data.frame(clean)

clean$GPS.coordinates.long <- as.numeric(as.character(clean$GPS.coordinates.long))
clean$GPS.coordinates.lat <- as.numeric(as.character(clean$GPS.coordinates.lat))

clean  <- clean[!is.na(clean$GPS.coordinates.long) & !is.na(clean$GPS.coordinates.long), ]


clean$number_of_mice <- as.numeric(as.character(clean$number_of_mice))
clean$number_of_traps <- as.numeric(as.character(clean$number_of_traps))


## assigning areas
clean$area <-
    ifelse(grepl("Emanuel", clean$collection_team),
           "F",
           ifelse(grepl("Irina", clean$collection_team),
                  "D",
                  ifelse(grepl("Jana", clean$collection_team),
                         "A",
                         ifelse(grepl("Gregor", clean$collection_team),
                                "C",
                                ifelse(grepl("Paula|Victor-Sascha-Naomi",
                                             clean$collection_team),
                                       "B", 
                                       ifelse(grepl("Sascha", clean$collection_team),
                                              "E", "No"))))))



clean.caught <- clean[clean$number_of_mice>0, ]

clean$Y_N<- ifelse(clean$number_of_mice>0, "Mm", "None")

clean$MperT <- clean$number_of_mice/clean$number_of_traps

write.csv(clean.caught, "Cleaned_HMHZ_2016_MmOnly.csv")

write.csv(clean, "Cleaned_HMHZ_2016_All.csv")


########### Overall #######

## how many mice
sum(clean$number_of_mice)

## how many trap nights
sum(clean$number_of_traps)

## trapping efficiency
sum(clean$number_of_mice)/sum(clean$number_of_traps)

## how many localities with mice
table(clean$number_of_mice>0)

########### the same divided by areas

## number of mice
num.mice.area <- tapply(clean$number_of_mice, clean$area, sum)

## trap nights
num.set.area <- cbind(tapply(clean$number_of_traps, clean$area, sum))

## localities with mice
loc.by.area <- tapply(clean$number_of_mice, clean$area,
                      function (x) table(x>0))

## trapping succes
eff.by.area <- by(clean, clean$area, function (x)
    sum(x$number_of_mice)/sum(x$number_of_traps))


team.table <- as.data.frame(cbind(num.mice.area, num.set.area,
                                  round(eff.by.area, 3),
                                  do.call(rbind, loc.by.area)))

colnames(team.table) <- c("mice", "traps", "trap_success",
                          "neg.local", "pos.local")

team.table$loc_success <- round(team.table$pos.local/
                                    (team.table$neg.local +
                                         team.table$pos.local), 3)


write.table(team.table, "team_table.csv", sep=";")


library(ggplot2)

## histogram of catches in localities 

png("figures/mice_per_locality.png", res=300, width=2400, height=1600)
ggplot(clean, aes(number_of_mice, fill=area)) +
    geom_histogram(binwidth=1) +
        scale_y_continuous("number of localites") +
            scale_x_continuous("number of mice caught") +
                ggtitle("distribution of trapping success at localities")
dev.off()

## number of traps
png("figures/number_traps_set.png", res=300, width=2400, height=1600)
ggplot(clean, aes(area, number_of_traps)) +
    geom_boxplot() +
        scale_y_continuous("number of traps") +
            scale_x_discrete("trapping area") +
                geom_jitter(width=0.1, color="red") +
                    ggtitle("distribution of traps set at localites")
dev.off()

ggplot(clean, aes(area, number_of_traps)) +
    geom_violin()+
        geom_jitter(width=0.1, color="red")


## number of mice
png("figures/number_mice_caught.png", res=300, width=2400, height=1600)
ggplot(clean, aes(area, number_of_mice)) +
    geom_boxplot()+
        scale_y_continuous("number of mice") +
            scale_x_discrete("trapping area") +
                geom_jitter(width=0.1, color="red") +
                    ggtitle("distribution of mice caught at localites")
dev.off()


## mouse per trap
png("figures/mice_per_trap.png", res=300, width=2400, height=1600)
ggplot(clean, aes(area, MperT)) +
    geom_boxplot() +
        scale_y_continuous("proportion of trapping success") +
            scale_x_discrete("trapping area") +
                geom_jitter(width=0.1, , color="red")+
                    ggtitle("trapping success mice cought per trap")
dev.off()


## mouse per trap
png("figures/mice_per_trap_log.png", res=300, width=2400, height=1600)
ggplot(clean, aes(area, MperT)) +
    geom_boxplot() +
        scale_y_log10("proportion of trapping success (log10)") +
            scale_x_discrete("trapping area") +
                geom_jitter(width=0.1, , color="red") +
                    annotation_logticks(sides="l")+
                        ggtitle("trapping success (mice caught per trap)")
dev.off()


ggplot(clean, aes(area, MperT)) +
    geom_violin()+
        geom_jitter(width=0.1, color="red") 

## number of traps vs trapping success per trap
ggplot(clean, aes(number_of_traps, MperT, color=area)) +
    geom_jitter() 


## trapping taktics differ per area?
ggplot(clean, aes(number_of_traps, MperT, color=area)) +
    geom_jitter() +
        stat_smooth(se=FALSE)


png("figures/mice_vs_traps.png", res=300, width=2400, height=1600)
ggplot(clean, aes(number_of_mice, number_of_traps, color=area))+
           scale_y_continuous("number of traps") +
               scale_x_continuous("number of mice") +
                   geom_jitter()+
                       stat_smooth(method="lm", se=FALSE) +
                           ggtitle("mice vs traps for all localities")
dev.off()

png("figures/mice_vs_traps_successLoc.png",
    res=300, width=2400, height=1600)
ggplot(subset(clean, number_of_mice>0),
       aes(number_of_mice, number_of_traps, color=area))+
           scale_y_continuous("number of traps") +
               scale_x_continuous("number of mice") +
                   geom_jitter()+
                       stat_smooth(method="lm", se=FALSE) +
                           ggtitle("mice vs traps for localities with mice only")
dev.off()


area <- get_map(location =
                    c(min(clean$GPS.coordinates.lat),
                      min(clean$GPS.coordinates.long),
                      max(clean$GPS.coordinates.lat),
                      max(clean$GPS.coordinates.long)),
                source = "stamen", maptype="toner-lite")

png("figures/Ludo.png", units = "in", width = 6, height = 7, res=300)
ggmap(area,  zoom = 15) +
    geom_point(data = subset(ludo, !is.na(ludo$HIX)),
               aes(Longitude, Latitude, color=HIX), size=2, alpha=0.5) +
                   scale_color_gradient("Hybrid\nindex", high="red",low="blue")
dev.off()
