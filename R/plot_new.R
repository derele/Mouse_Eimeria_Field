library(ggmap)
library(png)
library(grid)


raw <- read.csv("locations_table_corrected.csv")

LOC.catched <- raw[rowSums(raw[, c("female", "male", "juv")],
                           na.rm = TRUE)>0, ]

write.csv(LOC.catched, "/home/ele/Dropbox/field_sampling/catched_2014.csv")


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


LOC.catched$day_traps_removed <- as.character(
    as.Date(LOC.catched$day_traps_removed, "%d.%m.%y"),  "%m/%d/%y")

GEN <- read.csv("genotypes2014.csv")


mice <- merge(GEN, LOC.catched, by.x = "Code", by.y = "CODE_CZ_DB")

collapsed.GT <- apply(mice[, 7:14], 1, paste, collapse = "/")

hyb.d <- sapply(collapsed.GT, function (x){
    dom <- nchar(x) - nchar(gsub("d", "", x))
    mus <- nchar(x) - nchar(gsub("m", "", x))
    dom/(dom+ mus)
})


mice <- cbind(mice, hyb.d)

mice <- mice[!is.na(mice$hyb.d),]


library(plotrix)

area.pic <- readPNG("ggmapTemp.png")

clustered.coord <- cluster.overplot(mice$longitude, mice$latitude)# , away = c(0.02,0.02))

pdf("HZ_BR_2014_cluster_on_osm.pdf")
ggplot(mice, aes(clustered.coord$x, clustered.coord$y, color=hyb.d)) +
    annotation_custom(rasterGrob(area.pic, width=unit(1,"npc"), height=unit(1,"npc")),
                      -Inf, Inf, -Inf, Inf) +
    xlim(12.3, 14.2) +
    ylim(52.3, 53.5) +
    geom_point(size=1.5)+
    scale_color_gradient(low="red", high = "blue")
dev.off()


    
