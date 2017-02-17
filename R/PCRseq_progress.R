## Super map of PCR results 

## In the table PCR2016.csv, include the code corresponding to each ID_mouse (present in the file HZ16_Mice_18-0716_dissections.csv)

## Get codes associated to mice Id
codes.2016 <- read.csv("~/git_projects/Mouse_Eimeria_Databasing/raw_data/HZ16_Mice_18-07-16_dissections.csv")

## Copy the columns that I want from the original file and then eliminate the raws that contain information that won't be used
codes.2016 <- codes.2016[c(4,5)]
codes.2016 <- codes.2016[-c(1:5),]

## Add GPS coordinate for each code (locality)
gps.coord16 <- read.csv("~/git_projects/Mouse_Eimeria_Databasing/output_data/all_clean_localities.csv")

#Delete first column
gps.coord16 <- gps.coord16[-1]

## First update the file
pcr.16 <- read.csv("~/git_projects/Mouse_Eimeria_Databasing/raw_data/PCR_2016.csv") 

##merge the pcr results with the codes per individuals
pcr_codes_2016 <- merge(pcr.16, codes.2016, by= intersect(names(pcr.16),names(codes.2016)))

##merge the gps coords to pcr_codes per individuals
Total.matrix <- merge(pcr_codes_2016, gps.coord16, by= intersect(names(pcr_codes_2016),names(gps.coord16)), all=TRUE)

##Eliminate NA 
Total.matrix <- na.omit(Total.matrix)

## Include Hybrid Index to the table
##later...

## Plot in the map the PCR results 
##install.packages("ggmap", dependencies=TRUE, repos="http://cran.rstudio.com/")
#install.packages("ggmap", type = "source")
library(ggmap)
area <- get_map(location =
                  c(min(Total.matrix$Longitude-0.3),
                    min(Total.matrix$Latitude-0.3),
                    max(Total.matrix$Longitude+0.3),
                    max(Total.matrix$Latitude+0.3)),
                source = "stamen", maptype="toner-lite")

#png("Fakemap.png", units = "in", width = 8, height = 8, res=300)
library(ggrepel)

uniqueLoc16 <- unique(Total.matrix[c(1,10,11)])

ggmap(area)+
    geom_point(data = Total.matrix,  aes(Longitude, Latitude, color=as.factor(COCE_DNA)), size=5, alpha=0.6)+
    scale_color_manual(values=c("red", "darkgreen"), name="Status DNA", label=c("not extracted","extracted"))+
    geom_label_repel(data= uniqueLoc16, aes(Longitude, Latitude, label=Code),
                     arrow = arrow(length = unit(0.01, 'npc')),
                     force = 10)
