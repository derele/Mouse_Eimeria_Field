## November 2017, Alice Balard
library(ggplot2)
library(ggmap)
library(ggrepel)
library(data.table)

source("HMHZ_Functions.R")

MiceTable <- read.csv("../raw_data/MiceTable_2014to2017.csv")
TrapTable <- read.csv("../raw_data/TrapTable_2014to2017.csv")

## Remove the Poland data
MiceTable <- MiceTable[which(MiceTable$Longitude < 17),]

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
MiceTable$Latitude <- round(MiceTable$Latitude, 2)
MiceTable$Longitude <- round(MiceTable$Longitude, 2)
# ***************************************************************

# If the same cluster was sampled several time one year, we sum the mice caught
aggdata <- aggregate(x = TrapTable[c("Number_mus_caught", "Number_traps_set")],
                     by = TrapTable[c("Latitude", "Longitude", "Year")], FUN = sum)

names(aggdata)[names(aggdata) == "x"] <- "Number_mice_caught"

# add a counter for the localities:
aggdata <- as.data.table(aggdata)[, count := seq(.N), by = c("Latitude", "Longitude")][]

# density of mice in 100%
aggdata$density <- aggdata$Number_mus_caught / aggdata$Number_traps_set *100

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

# density of mice added to MiceTable
TotalTable <- merge(MiceTable, aggdata, all = TRUE)

TotalTable <- TotalTable[TotalTable$Longitude < 17,]
# ***************************************************************

# Sex - localities 
ggplot(TotalTable[TotalTable$Sex %in% c("M", "F"), ], aes(x = Longitude, y = Latitude, fill = Sex)) +
  theme_classic(base_size = 18) +
  theme(axis.title.y=element_blank(), axis.title.x=element_blank()) +
  ggtitle(label = "Cluster of males and female per localities") +
  theme(legend.title=element_blank()) +
  scale_fill_manual(values = c("pink", "blue")) +
  geom_density2d(aes(color = Sex), size = 1) +
  geom_point(col = "black", size = 4, pch = 21)
  
# how many localities successful per year in 2016 and 2017? 
ggplot(aggdata[aggdata$Year %in% c(2016,2017),], aes(x=factor(Year)))+
  geom_bar(aes(fill = is.success), stat="count", width=0.7, color="black")+
  theme_classic(base_size = 18) +
  theme(axis.title.y=element_blank(), axis.title.x=element_blank()) +
  ggtitle(label = "Number of localities trapped") +
  theme(legend.title=element_blank()) +
  scale_fill_manual(values = c("darkgrey", "green"))

# Host density - sex
ggplot(TotalTable[TotalTable$Sex %in% c("M", "F"),], aes(x = density))+
  geom_bar(aes(fill = Sex), stat = "count", binwidth = 3, col = "white") +
  theme_classic()

# Host density - hybrid index
ggplot(TotalTable[!is.na(TotalTable$HI),], aes(x = HI, y = density))+
  geom_point(pch = 21)+
  geom_smooth(col = "red") +
  theme_classic()

# how many mice where caught every year?
aggmice <- aggregate(x = TrapTable[c("Number_mus_caught")],
                     by = TrapTable[c("Year")], FUN = sum)

ggplot(aggmice, aes(x = factor(Year), y = Number_mus_caught))+
  geom_bar(stat="identity", width=0.7, color="black", fill = "orange")+
  theme_classic(base_size = 18) +
  theme(axis.title.y=element_blank(), axis.title.x=element_blank()) +
  ggtitle(label = "Number of mice caught per year") +
  geom_text(aes(label=Number_mus_caught), vjust=1.6, color="black", size=10)+
  theme(legend.title=element_blank()) 

# Number of localities with 1, 2, 3 or more mice caught over the years
ggplot(aggdata, aes(x = Number_mus_caught)) +
  geom_histogram(fill="black", col="grey", binwidth = 1) +
  ggtitle(label = "Histogram of # mice caught per locality")+
  theme_classic(base_size = 18) +
  labs(x = "Number of mice caught")

# all the possible worms: (to improve with Jenny)
TotalTable$worms_count <- rowSums(TotalTable[c("Unknown_cecum","Unknown_colon",	"Unknown_SI",	"Syphacia_cecum_colon",	"Cysticercus_tenia_teniformis_liver",	
                                             "Catotenia_pusilla_SI",  "Mastophorus_stomach_SI",	"Heterakis_spumosa_colon_cecum",	"Tapeworm_SI",	
                                             "Trichuris_cecum",	"Aspiculuris_cecum",  "Aspiculuris_colon",	"Hymenolepis_SI",	"Rodentolepis_liver_digtract",
                                             "Mesocoides_body_cavities",	"Mesocoides_lungs",  "Heigmosomoides_polyguis",	"H_diminita.")],
                                 na.rm = TRUE)
TotalTable$worms_presence <- 0
TotalTable$worms_presence[TotalTable$worms_count !=0] <- 1

aggworms <- aggregate(x = TotalTable[c("worms_presence")],
                      by = TotalTable[c("Latitude", "Longitude", "Year")], FUN = sum)

aggdata <- merge(aggdata, aggworms, all = TRUE)
aggdata$worms_prevalence <- aggdata$worms_presence / aggdata$Number_mus_caught *100

ggplot(aggdata, aes(x = density, y = worms_prevalence)) +
  geom_smooth(fill = "orange", color = "brown") +
  geom_jitter(shape = 21, fill = "brown", size = 3) +
  coord_cartesian(ylim = c(0,99)) +
  theme_classic(base_size = 18) +
  labs(x = "Mice density")

# prevalence at a locality depending on the density

#***************************

# Map
margin <- 0.5
area <- get_map(location =
                  c(min(na.omit(MiceTable$Longitude) - margin),
                    min(na.omit(MiceTable$Latitude) - margin),
                    max(na.omit(MiceTable$Longitude) + margin),
                    max(na.omit(MiceTable$Latitude) + margin)),
                source = "stamen", maptype="toner-lite",
                zoom = 7)

#plot the map :
ggmap(area) +
  geom_point(data = MiceTable, shape = 21, size = 4,
             aes(Longitude, Latitude, fill = as.factor(Year))) +
  theme(legend.text=element_text(size=20)) +
  guides(fill=guide_legend(title=""))


#***************************
## To finish

# Calculate HI per location with ratio of alleles

# link eimeria and prevalence


# weight --> age

summary(MiceTable$Body_weight)


MiceTable$subspecies[MiceTable$HI < 0.1] <- "Mmd"
MiceTable$subspecies[MiceTable$HI >= 0.9] <- "Mmm"
#MiceTable$subspecies[MiceTable$HI ] <- "Mmm"


ggplot(MiceTable, aes(x = Body_weight, y = Body_length, fill = subspecies, order = order)) +
  geom_point(col = "black", pch = 21, size = 5) +
  geom_point(data = MiceTable[MiceTable$subspecies %in% c("Mmm", "Mmd"), ], col = "black", pch = 21, size = 5) +
  scale_fill_manual(values = c("blue", "red", "grey")) +
  theme_classic(base_size = 18) +
  facet_grid(. ~ Sex)

ggplot(MiceTable[MiceTable$subspecies %in% c("Mmm", "Mmd"), ], aes(x = Body_weight, y = Body_length, fill = subspecies, order = order)) +
  geom_density2d(aes(color = subspecies), size = 2) +
  scale_fill_manual(values = c("blue", "red", "grey")) +
  theme_classic(base_size = 18) 

# body condition index to correlate with parasite data correlate with hybrid index

