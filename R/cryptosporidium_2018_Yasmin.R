## First assessment of cryptosporidium in HMHZ
## 2018 BSc Yasmin Schickel
library(ggplot2)
library(ggmap)
library(viridis) # pretty colors
library(ggrepel) # to make labels on plot

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
                  c(min(Yasdata$Longitude - 0.5),
                    min(Yasdata$Latitude - 0.5),
                    max(Yasdata$Longitude + 0.5),
                    max(Yasdata$Latitude + 0.5)),
                source = "stamen", maptype= "toner-lite",
                zoom = 8)

#plot the map :

# little trick for ggrepel
rownames(fullYasminData) <- fullYasminData$Mouse_ID


ggmap(area) +
  geom_point(data = Yasdata, 
             aes(x = Longitude, y = Latitude, 
                 fill = log10(Nbr_oocyst_MEAN + 1)),
             pch = 21, size = 4, alpha = 0.8) +
  scale_fill_viridis(option="inferno") +
  geom_label_repel(data = Yasdata, 
                   arrow = arrow(length = unit(0.01, "npc"), 
                                 type = "closed", ends = "last"),
                   force = 10,
                   aes(x = Longitude, y = Latitude, 
                       label = rownames(Yasdata)),
                       size = 3)
