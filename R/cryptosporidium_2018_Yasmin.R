## First assessment of cryptosporidium in HMHZ
## 2018 BSc Yasmin Schickel
library(ggplot2)
library(ggmap)
library(viridis) # pretty colors
library(ggrepel) # to make labels on plot
library(png) # to add image on plot (mice)

# Import data
yasminData <- read.csv("../data/Cryptosporidium/cryptoDetection2018.csv") 
mice2018data <- read.csv("../data/Field_data/HZ18_September_Mice_Dissection.csv")

# Homogenize the column Mouse_ID
yasminData$Mouse_ID <- gsub(pattern = "IL ", replacement = "AA_0", x = yasminData$Sample_ID) # gsub change a pattern for another

# Add information we got from the field by merging
fullYasminData <- merge(yasminData, mice2018data)

# Let's make a map
Yasdata <- fullYasminData[fullYasminData$qPCR_results == "TRUE",]

# get a map
area <- get_map(location =
                  c(min(Yasdata$Longitude - 1),
                    min(Yasdata$Latitude - 1),
                    max(Yasdata$Longitude + 1),
                    max(Yasdata$Latitude + 1)),
                source = "stamen", maptype= "toner-lite",
                zoom = 8)

#plot the map :

# little trick for ggrepel
rownames(Yasdata) <- Yasdata$Mouse_ID

ggmap(area) +
  scale_fill_viridis(option="inferno") +
  geom_label_repel(data = Yasdata,min.segment.length = .7,
                   arrow = arrow(length = unit(0.01, "npc"), 
                                 type = "closed", ends = "last"),
                   force = 10,
                   aes(x = Longitude, y = Latitude, 
                       label = rownames(Yasdata)),
                       size = 3) +
  geom_point(data = Yasdata, 
             aes(x = Longitude, y = Latitude, 
                 fill = log10(Nbr_oocyst_MEAN + 1)),
             pch = 21, size = 4) 

## The most pretty plot EVER

g <- ggmap(area) +
  scale_fill_viridis(option="inferno") +
  geom_label_repel(
    aes(x = Longitude, y = Latitude, 
        fill = Nbr_oocyst_MEAN,
        label = rownames(subset(Yasdata, Longitude < 14))),
    data          = subset(Yasdata, Longitude < 14),
    color         = "grey",
    nudge_x       = 3 - subset(Yasdata, Longitude < 14)$Longitude,
    segment.size  = 0.2,
    direction     = "y",
    hjust         = 0,
    segment.colour = "black"
    ) +
  geom_label_repel(
    aes(x = Longitude, y = Latitude, 
        fill = Nbr_oocyst_MEAN,
        label = rownames(subset(Yasdata, Longitude >= 14))),
    data          = subset(Yasdata, Longitude >= 14),
    color         = "lightgrey",
    nudge_x       = 3 + subset(Yasdata, Longitude >= 14)$Longitude,
    segment.size  = 0.2,
    direction     = "y",
    hjust         = 0,
    segment.colour = "black"
  ) 

mice <- readPNG("mousepic.png")

for (i in 1:nrow(Yasdata)) { 
  g <- g + inset_raster(mice, 
                        xmin = Yasdata$Longitude[i]- 0.18, 
                        xmax = Yasdata$Longitude[i]+ 0.18,
                        ymin = Yasdata$Latitude[i]- 0.05, 
                        ymax = Yasdata$Latitude[i]+ 0.05, 
                        interpolate = TRUE) }
g
