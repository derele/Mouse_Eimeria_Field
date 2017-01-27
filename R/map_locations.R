if(!exists("loc")){
    source("R/locality_import.R")
}
    
library(ggmap)

#######################################
area.B <- get_map(location= c(min(loc$Longitude)-0.3,
                              min(loc$Latitude)-0.3,
                              max(loc$Longitude)+0.3,
                              max(loc$Latitude)+0.3),
                  source = "stamen", maptype="toner-lite")

##                        min(loc$Longitude),
##                        max(loc$Latitude),

pdf("figures/Sanity2016.pdf", width=50, height=50)
ggmap(area.B) +
        geom_text(data = loc,
                   aes(Longitude, Latitude, label=Code)) 
dev.off()


