# Clean qPCR data from start
# Eimeria gene / Mouse gene, triplicates

# preparation : 
# convert all xls to csv files: for i   in *.xlsx; do  libreoffice --headless --convert-to csv "$i" ; done
# stack them together: cat *.csv >> allBavaria2015.csv (done in bash). 
# remove manually the useless lines

library(dplyr)
library(ggplot2)

# Import data:
bav2015 <- read.csv("allBavaria2015.csv")

# NB only CEWE!!!

# how many mice?
length(unique(bav2015$Name)) 

# remove those more than triple
table(bav2015$Name)[table(bav2015$Name) > 6]

# CEWE_SK_2973 CEWE_SK_2977 
# 12           12 
# CEWE_SK_3024 CEWE_SK_3027 
# 12           12 
# CEWE_SK_3032 CEWE_SK_3033 
# 18           12 
# CEWE_SK_3039 CEWE_SK_3047 
# 18           12 
# CEWE_SK_3048 CEWE_SK_3050 
# 12           12 
# CEWE_SK_3062 CEWE_SK_3063 
# 12           12 
# CEWE_SK_3064 CEWE_SK_3065 
# 12           12 
# CEWE_SK_3068 
# 12 

head(bav2015)

# format data to have "tissue" and "Mouse_ID" column
bav2015$Name

x <- strsplit(as.character(bav2015$Name), "_", 1) # usefull function here :)
x
bav2015$tissue <- sapply(x, "[", 1) # means "take the first (1) element of list ("[") in x (x)
bav2015$Mouse_ID <- paste(sapply( x, "[", 2), sapply( x, "[", 3), sep = "_")
rm(x)

## Function to calculate deltaCtMminusE
calculateDeltaCt <- function(df, mergeBy){
  sumDataMouse <- df[df$Target.SYBR %in% "mouse",]
  sumDataEimeria <- df[df$Target.SYBR %in% "eimeria",]
  mergedData <- merge(sumDataEimeria, sumDataMouse,
                      by = mergeBy)
  mergedData <- mergedData[!duplicated(mergedData[,c('Mouse_ID','tissue')]),]
  mergedData$deltaCt_MminusE <- as.numeric(as.character(mergedData$Ct.Mean.SYBR.y)) -
    as.numeric(as.character(mergedData$Ct.Mean.SYBR.x)) # DeltaCt MOUSE minus EIMERIA
  return(mergedData)
}
  
bav2015 <- calculateDeltaCt(bav2015, mergeBy = c("Mouse_ID", "tissue"))

head(bav2015)

# So far we work JUST on CEWE, NB!

# Some plot to see

ggplot(bav2015, aes(x = Mouse_ID, y = deltaCt_MminusE)) +
  geom_point(aes(col = deltaCt_MminusE > -5))

# add HI :) but we have just 13 positive so for model would be too few...
miceTable <- read.csv("../../../MiceTable_fullEimeriaInfos_2014to2017.csv")

fullBav2015 <- merge(bav2015, miceTable, by = "Mouse_ID")
# and names match all so we are SUPER HAPPY

# for plot...
mypngfile <- download.file('http://clipart-library.com/images/piqKdL7i9.png', destfile = 'mypng.png', mode = 'wb')
library(png)
mypng <- readPNG('mypng.png')

# plot along HI
ggplot(fullBav2015[fullBav2015$deltaCt_MminusE > -5, ], 
       aes(x = HI, y = deltaCt_MminusE)) +
  geom_point() +
  geom_smooth() + 
  theme_classic() +
  annotate("text", x = 0.1, y = -1, 
           label = "this is just CEWE for Bavaria 2015 \n and not enough points to make any conclusion...") +
  annotation_raster(mypng, ymin = -2.3,ymax= -1.5,xmin = 0.04,xmax = 0.1)

