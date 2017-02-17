## Map of Enas' infection experiment
source("/home/alice/git_projects/Mouse_Eimeria_Databasing/R/Concatenated_data.R")

MiceInfID <- c("AA_0064", "AA_0068", "AA_0070", "AA_0139", "AA_0088")

formap <- dissection.2016[dissection.2016$ID_mouse %in% MiceInfID,c(4,5)]
formap <- merge(All_loc,formap)

## one general map
areainf <- get_map(location =
                    c(min(formap$Longitude)-0.5,
                      min(formap$Latitude)-0.5,
                      max(formap$Longitude)+0.5,
                      max(formap$Latitude)+0.5),
                source = "stamen", maptype="toner-lite")


## install.packages("ggmap", type = "source")
#png("/home/alice/git_projects/Mouse_Eimeria_Databasing/output_data/Concatenated_map_feb2017.png", units = "in", width = 6,
 #   height = 7, res=300)
ggmap(areainf)+
       geom_jitter(data = formap,
               aes(Longitude, Latitude), size=4, alpha=1, color = "green")
 #   scale_color_gradient("Hybrid\nindex", high="red",low="blue")
#dev.off()                         
