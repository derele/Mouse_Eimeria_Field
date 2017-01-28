library(ggmap)
library(png)
library(grid)


raw <- read.csv("Brandenburg_sampling_2014.csv")
raw <- raw[, 4:ncol(raw)-2]


LOC <- raw[, c("latitude", "longitude", "female", "male", "juveniles")]

LOC.catched <- LOC[rowSums(LOC[, c("female", "male", "juveniles")],
                           na.rm = TRUE)>0, ]

LOC.catched$juveniles[is.na(LOC.catched$juveniles)] <- 0

LOC.catched$mice <- rowSums(LOC.catched[3:5])

LOC.catched$longitude <- gsub(",", ".", LOC.catched$longitude)

LOC.catched$latitude <- gsub(",", ".", LOC.catched$latitude)

write.csv(LOC.catched, "/home/ele/Dropbox/field_sampling/catched_2014.csv")


LOC.catched$latitude <- as.numeric(LOC.catched$latitude)
LOC.catched$longitude <- as.numeric(LOC.catched$longitude)

LOC.catched$mice <- as.integer(LOC.catched$mice)


max(LOC.catched$longitude)
min(LOC.catched$longitude)
max(LOC.catched$latitude)
min(LOC.catched$latitude)

area <- get_map(location = c(12.3, 52.3, 14.2, 53.5),
                source = 'osm')
area

area.pic <- readPNG("ggmapTemp.png")


png("area_dots.png", width = 20, height = 16, units = 'in', res = 300)
ggplot(LOC.catched, aes(LOC.catched$longitude, LOC.catched$latitude, size = LOC.catched$mice)) +
    annotation_custom(rasterGrob(area.pic, width=unit(1,"npc"), height=unit(1,"npc")),
                      -Inf, Inf, -Inf, Inf) +
    xlim(12.3, 14.2) +
    ylim(52.3, 53.5) +
    geom_point(color="black", size=LOC.catched$mice + 3)+ 
    geom_point(color="red", fill="darkred", size=LOC.catched$mice+2)+
    scale_size("number of mice")

dev.off()


locations <- raw[, c("location", "address", "location_code", "latitude", "longitude", "day_traps_removed", "female", "male", "juveniles")]

locations <- locations[rowSums(LOC[, c("female", "male", "juveniles")],
                               na.rm = TRUE)>0, ]

locations$longitude <- gsub(",", ".", locations$longitude)
locations$latitude <- gsub(",", ".", locations$latitude)

locations[is.na(locations)] <- 0

locations$day_traps_removed <- as.character(
    as.Date(locations$day_traps_removed, "%d.%m.%y"),  "%m/%d/%y")

rownames(locations) <- NULL

total_mice_loc <- rowSums(locations[, c("female", "male", "juveniles")])

mice <- locations[rep(seq_len(nrow(locations)), times=total_mice_loc),]

rownames(mice) <- NULL

locations$location_code <- as.character(locations$location_code)

sex.vec <- unlist(by(locations, locations$location_code, function (x){
    c(rep("M", times = x$male),
      rep("F", times = x$female),
      rep("JUV", times = x$juveniles))
}))


mice <- mice[order(mice$location_code),]

mice$sex <- sex.vec

mice <- mice[, !names(mice)%in%c("female", "male", "juveniles")]

write.csv(locations, "locations_table.csv", row.names=FALSE)

write.csv(mice, "mice_table.csv", row.names=FALSE)
