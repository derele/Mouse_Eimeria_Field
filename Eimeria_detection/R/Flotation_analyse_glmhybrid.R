#########################
library(ggplot2)
library(gridExtra)
library(reshape2)
library(scales)
library(RCurl)
library(ggrepel)
library(ggmap)
library(devtools)
library(MASS)

#########################
## Import :
data2015 <- read.csv("../data_clean/Results_flotation_and_PCR_2015_CLEAN_temp.csv")
data2016 <- read.csv("../data_clean/Results_flotation_and_PCR_2016_CLEAN.csv")

data <- rbind(data2015, data2016)

#########################
## Add the more recent HI :
source("https://raw.githubusercontent.com/alicebalard/EimeriaHMHZ_various/master/HMHZ_functions.R")
datafull <- merge(data, getGenDF(), by = "Mouse_ID")

## NB later! Add haplotype (A or B)

## How many mice caught :
length(unique(datafull$Mouse_ID))

## Which are the duplicated samples?
datafull[which(duplicated(datafull$Mouse_ID)),]$Mouse_ID

## How many sampling sites :
length(unique(datafull$Code))

#########################
## Calculate prevalence by sampling site :
data.agg <- aggregate(na.omit(datafull)$OPG, by = list(na.omit(datafull)$Code), 
                      FUN = sum, na.rm=TRUE)

length(data.agg[which(data.agg$x != 0), ]$Group.1) / length(data.agg$Group.1) *100

#########################
## Number of sites sampled :
datafull$Code <- factor(datafull$Code)
length(levels(datafull$Code))
datafull$HI <- as.numeric(datafull$HI)

sample_data <- function(datafull){
  ## Select 1 individual mouse per farm randomly :
  # a. iterate over the unique values of "Code"
  # b. find the row indices for each value
  # c. select one row index at random 
  ind <- sapply(unique(datafull$Code) , function(x) sample(which(datafull$Code==x) , 1))
  dataRan <- datafull[ind, ]
  dataRan$HI <- as.numeric(dataRan$HI)
  dataRan$OPG <- as.numeric(dataRan$OPG)
  
  ## Cheat for glm.hyrid (needs to be inplemented!!)
  #dataRan[nrow(dataRan)+1,]$HI <- 1 
  #dataRan[nrow(dataRan),]$OPG <- 0
  #is.factor(dataRan$Code)
  
  #########################
  ## Run glm.hybrid on the data :
  # install_github("alicebalard/Parasite_Load")
  
  # glm.hybrid::glm.hybrid(formula = OPG ~ HI * Code, data = dataRan, alpha.along = "HI")
  # Fail!!!!!
  
  ## Compare :
  model <- glm.nb(formula = OPG ~ HI, data = dataRan)
  summary(model)
  
  ## Quick plot :
  plot(x = dataRan$HI, y = dataRan$OPG)
  hi.along <- seq(0, 1, 0.1)
  lines(hi.along, predict(model, data.frame(HI = hi.along),
                          type = "response"))
  
}

#########################
## Bootstrap this 100 times, summarize :
sample_data(datafull)

#########################
## Map all samples :
datafull$HI <- as.numeric(datafull$HI)

# get a map
margin <- 2

area <- get_map(location =
                  c(min(datafull$Longitude - margin),
                    min(datafull$Latitude - margin),
                    max(datafull$Longitude + margin),
                    max(datafull$Latitude + margin)),
                source = "stamen", maptype="toner-lite",
                zoom = 7)

#plot the map :
map <- ggmap(area) +
  geom_point(data = datafull, shape = 21,
             aes(Longitude, Latitude, fill = HI, size = OPG), alpha = 0.5) + # set up the points
  scale_fill_gradient("Hybrid\nindex", high="red",low="blue") + # set up the HI colors
  scale_size_continuous(range = c(2, 10), name = "Oocyst per gram")
map

# Export :
pdf("../figures/Map2015-2016_temp.pdf")
plot(map)
dev.off()
