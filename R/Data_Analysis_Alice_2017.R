## November 2017, Alice Balard
library(ggplot2)
library(ggmap)
library(ggrepel)
library(data.table)

source("HMHZ_Functions.R")

MiceTable <- read.csv("../raw_data/MiceTable_2014to2017.csv")
TrapTable <- read.csv("../raw_data/TrapTable_2014to2017.csv")

## Let's split the localities VISITED and the localities where MICE WERE TRAPPED (before 2016)

# ***************************************************************
## Cluster localities : how many times each LOCALITY was sampled?

# pairsloc <- pairwise.cluster.loc(trapsTOT) # errors with this model, to discuss
# loctime <- data.frame(table(apply(pairsloc, 2, sum)))
# loctime$Var1 <- as.numeric(as.character(loctime$Var1)) + 1 
# names(loctime) <- c("Repeatition over the years", "# localities")

measure(lon1 = 13.04390, lat1 = 51.93000, lon2 = 13.64000, lat2 = 52.71690) - measure(lon1 = 13.04, lat1 = 51.93, lon2 = 13.64, lat2 = 52.72)
measure(lon1 = 13.6761, lat1 = 52.4978, lon2 = 13.68, lat2 = 52.50)

# rounded to about 700 meters
TrapTable$Longitude <- round(TrapTable$Longitude, 2)
TrapTable$Latitude <- round(TrapTable$Latitude, 2)
# ***************************************************************

# If the same cluster was sampled several time one year, we sum the mice caught
aggdata <- aggregate(x = TrapTable[c("Number_mus_caught", "Number_traps_set")],
                     by = TrapTable[c("Latitude", "Longitude", "Year")], FUN = sum)

names(aggdata)[names(aggdata) == "x"] <- "Number_mice_caught"

# add a counter for the localities:
aggdata <- as.data.table(aggdata)[, count := seq(.N), by = c("Latitude", "Longitude")][]

# so you define new and old localities:
aggdata$is.new.loc <- "new"
aggdata[aggdata$count != 1, "is.new.loc"] <- "old"

# we define also successful and unsuccessful localities:
aggdata$is.success <- "failure"
aggdata[aggdata$Number_mus_caught != 0, "is.success"] <- "success"

# NB : in 2016 and 2017 we have the trapping attemps, not the other years

# how many localities (new and old) sampled per year? 
ggplot(aggdata[aggdata$Number_mus_caught != 0,], aes(x=factor(Year)))+
  geom_bar(aes(fill = is.new.loc), stat="count", width=0.7, color="black")+
  theme_classic(base_size = 18) +
  theme(axis.title.y=element_blank(), axis.title.x=element_blank()) +
  ggtitle(label = "Number of localities with mice where caught per year") +
  theme(legend.title=element_blank()) +
  scale_fill_manual(values = c("darkgreen", "gray"))

# hum I should have around 80 in 2015 according to Jenny's thesis...

# ***************************************************************

# how many localities successful per year in 2016 and 2017? 
ggplot(aggdata[aggdata$Year %in% c(2016,2017),], aes(x=factor(Year)))+
  geom_bar(aes(fill = is.success), stat="count", width=0.7, color="black")+
  theme_classic(base_size = 18) +
  theme(axis.title.y=element_blank(), axis.title.x=element_blank()) +
  ggtitle(label = "Number of localities with successful trapping per year") +
  theme(legend.title=element_blank()) +
  scale_fill_manual(values = c("darkgrey", "green"))
# ***************************************************************

# how many mice where caught every year?

aggmice <- aggregate(x = TrapTable[c("Number_mus_caught")],
                     by = TrapTable[c("Year")], FUN = sum)

ggplot(aggmice, aes(x=factor(Year), y = Number_mus_caught))+
  geom_bar(stat="identity", width=0.7, color="black", fill = "orange")+
  theme_classic(base_size = 18) +
  theme(axis.title.y=element_blank(), axis.title.x=element_blank()) +
  ggtitle(label = "Number of mice caught per year") +
  geom_text(aes(label=Number_mus_caught), vjust=1.6, color="black", size=10)+
  theme(legend.title=element_blank()) 
#***************************

# Number of localities with 1, 2, 3 or more mice caught over the years

ggplot(aggdata, aes(x = Number_mus_caught)) +
  geom_histogram(fill="black", col="grey", binwidth = 1) +
  ggtitle(label = "Histogram of # mice caught per locality")+
  theme_classic(base_size = 18) +
  labs(x = "Number of mice caught")

#***************************
# prevalence of worms depending on the number of mice caught

# all the possible worms: (to improve with Jenny)
MiceTable$worms_count <- rowSums(MiceTable[c("Unknown_cecum","Unknown_colon",	"Unknown_SI",	"Syphacia_cecum_colon",	"Cysticercus_tenia_teniformis_liver",	
                                             "Catotenia_pusilla_SI",  "Mastophorus_stomach_SI",	"Heterakis_spumosa_colon_cecum",	"Tapeworm_SI",	
                                             "Trichuris_cecum",	"Aspiculuris_cecum",  "Aspiculuris_colon",	"Hymenolepis_SI",	"Rodentolepis_liver_digtract",
                                             "Mesocoides_body_cavities",	"Mesocoides_lungs",  "Heigmosomoides_polyguis",	"H_diminita.")],
                                 na.rm = TRUE)

Worms <- MiceTable[!is.na(MiceTable$worms_count),c("Latitude", "Longitude", "Year", "worms_count")]
Worms$mice_caught <- 1

Worms$Latitude <- round(Worms$Latitude,2); Worms$Longitude <- round(Worms$Longitude,2)

Worms <- aggregate(x = Worms[c("worms_count", "mice_caught")],
          by = Worms[c("Latitude", "Longitude", "Year")], FUN = sum)

Worms$prevalence <- 0
Worms[Worms$worms_count != 0,]$prevalence <- 1

ggplot(Worms, aes(x = mice_caught, y = prevalence)) +
  geom_point(shape = 21, alpha = 0.1, fill = "brown", size = 10) +
  geom_smooth(fill = "orange", color = "brown") +
  coord_cartesian(xlim=c(0,20)) +
  ggtitle(label = "Prevalence of worms infection depending on the mice density",
          subtitle = "loess smoothing + 95%CI")+
  theme_classic(base_size = 18) +
  labs(x = "Mice prevalence")

#***************************

# Map
margin <- 0.5
area <- get_map(location =
                  c(min(na.omit(MiceTable$Longitude) - margin),
                    min(na.omit(MiceTable$Latitude) - margin),
                    max(na.omit(MiceTable$Longitude) + margin),
                    max(na.omit(MiceTable$Latitude) + margin)),
                source = "stamen", maptype="toner-lite",
                zoom = 6)

#plot the map :
ggmap(area) +
  geom_point(data = MiceTable, shape = 21, size = 3,
             aes(Longitude, Latitude, fill = as.factor(Year))) +
  theme(legend.text=element_text(size=20)) +
  guides(fill=guide_legend(title=""))

## To ask emanuel :
ggmap(area) +
  geom_point(data = MiceTable, shape = 21, size = 3,
             aes(Longitude, Latitude, fill = as.factor(HI)),show.legend=F) +
  theme(legend.text=element_text(size=20)) +
  guides(fill=guide_legend(position = NULL))
  
#***************************
## To finish

# Calculate HI per location with ratio of alleles

# link eimeria and prevalence