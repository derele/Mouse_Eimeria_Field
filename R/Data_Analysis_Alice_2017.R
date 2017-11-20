## November 2017, Alice Balard
library(ggplot2)
library(ggmap)
library(data.table)

source("HMHZ_Functions.R")

MiceTable <- read.csv("../raw_data/MiceTable_2014to2017.csv")
TrapTable <- read.csv("../raw_data/TrapTable_2014to2017.csv")

## Remove the Poland data
MiceTable <- MiceTable[which(MiceTable$Longitude < 17),]

## Remove the non Mus musculus (until further notice)
MiceTable$Species[is.na(MiceTable$Species)] <- "Mus musculus"

MiceTable <- MiceTable[which(MiceTable$Species == "Mus musculus"),]

# Check calculation of ALL HI
markers <- c("mtBamH", "Zfy2", "SRY1", "Y", "X332", "X347", "X65", "Tsx", "Btk",
             "Syap1", "Es1C", "Gpd1C", "Idh1C", "MpiC", "NpC", "Sod1C")

MiceTable <- get.HI.full(df = MiceTable, markers.col = markers)

# to bear in mind : jarda seems to have calculated HI on different sets of markers
table(MiceTable$HI == MiceTable$HI.calculated)

# HI per locality
agglocHI <- aggregate(x = MiceTable[c("collapsed.GT")],
          by = MiceTable[c("Latitude", "Longitude")], FUN = paste, sep = "/" )

agglocHI$HI <- get.HI(agglocHI$collapsed.GT)

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

# Map
HI.map(agglocHI, size = 4, alpha = 1, margin = 0.2,zoom = 7)
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

# NB : in 2016 and 2017 we have the trapping attemps, not the other years
# ***************************************************************

# density of mice added to MiceTable
TotalTable <- merge(MiceTable, aggdata, all = TRUE)

TotalTable <- TotalTable[TotalTable$Longitude < 17,]
# ***************************************************************

# body mass index as log body mass/log body length (Hayes et al. 2014)
TotalTable$BMI <- log(TotalTable$Body_weight) / log(TotalTable$Body_length)

# ***************************************************************

# how many localities (new and old) sampled per year? 

# barplot 
ggplot(data=data.frame(table(aggdata$Year, aggdata$is.new.loc)),
       aes(x=Var1, y=Freq, fill=Var2)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=Freq, label=Freq), vjust=1.6, 
            color="white", size=6)+
  scale_fill_brewer(palette="Paired")+
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

# ***************************************************************
# test potential linear correlation
library(devtools)
install_github("hadley/ggplot2")
install_github("ggobi/ggally")
library(GGally)

ggpairs(
  TotalTable, which(names(TotalTable) %in% c("Year","HI", "density", "BMI")),
  upper = list(
    continuous = wrap('cor', method = "spearman")
  ), 
  lower = list(
    continuous = 'cor'
  )
)

ggcorr(TotalTable[(names(TotalTable) %in% c("Year", "HI", "density", "BMI"))], 
       palette = "RdBu", label = TRUE, method = c("pairwise", "spearman")) +
  ggtitle("Pairwise correlations with Spearman coefficient")

# Sex - localities 
ggplot(TotalTable[TotalTable$Sex %in% c("M", "F"), ], aes(x = Longitude, y = Latitude, fill = Sex)) +
  theme_classic(base_size = 18) +
  theme(axis.title.y=element_blank(), axis.title.x=element_blank()) +
  ggtitle(label = "Cluster of males and female per localities") +
  theme(legend.title=element_blank()) +
  scale_fill_manual(values = c("pink", "blue")) +
  geom_density2d(aes(color = Sex), size = 1) +
  geom_point(col = "black", size = 4, pch = 21)

# Host density - sex
ggplot(TotalTable[TotalTable$Sex %in% c("M", "F"),], aes(x = density))+
  geom_bar(aes(fill = Sex), stat = "count", binwidth = 3, col = "white") +
  theme_classic()

# BMI - subspecies
TotalTable$subspecies[TotalTable$HI < 0.1] <- "Mmd"
TotalTable$subspecies[TotalTable$HI >= 0.9] <- "Mmm"

# Density plots 
ggplot(TotalTable[!is.na(TotalTable$subspecies),], aes(x = BMI, fill = subspecies)) +
  geom_density(alpha=.3) +
  scale_fill_manual(values = c("blue", "red"))

ggplot(TotalTable[!is.na(TotalTable$subspecies),], aes(x=subspecies, y = BMI, fill=subspecies)) +
  geom_boxplot() +
  guides(fill=FALSE) +
  scale_fill_manual(values = c("blue", "red")) +
  theme_classic()

# Host density - hybrid index
ggplot(TotalTable[!is.na(TotalTable$HI),], aes(x = HI, y = density))+
  geom_point(pch = 21)+
  geom_smooth(col = "red") +
  theme_classic()

ggplot(TotalTable[!is.na(TotalTable$HI),], aes(x = HI, y = density))+
  geom_point(pch = 21)+
  geom_smooth(col = "red", method = "lm") +
  theme_classic()

# ***************************************************************
# Parasites
WormsDF <- TotalTable[c("Mouse_ID", "Cysticercus", "Trichuris_muris", "Aspiculuris_tetraptera", "Syphacia_obvelata",
                        "Mastophorus_muris", "Heterakis_spumosa", "Mesocestoides", 
                        "Catenotaenia_pusilla", "Hymenolepis", "Oxyurids", "Mix_Syphacia_Aspiculuris",
                        "Heligmosomoides_polygurus", "Latitude", "Longitude")]

# ("Hymenolepis_microstoma", "Hymenolepis_diminiuta") to check, contains "TRUE" and "FALSE"

WormsDF <- na.omit(melt(WormsDF, id = c("Mouse_ID", "Longitude", "Latitude")))

aggworms <- aggregate(x = WormsDF["value"],
                      by = WormsDF["variable"], FUN = length)

ggplot(data=aggworms, aes(x=variable, y=value)) +
  geom_bar(stat="identity", fill = "#B983FF", col = "grey") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_text(aes(y=value, label=value), vjust=1.6, 
            color="black", size=6)

# prevalence at a locality depending on the density (only 2016, 2017)

WormsDF2 <- TotalTable[c("Cysticercus", "Trichuris_muris", "Aspiculuris_tetraptera", "Syphacia_obvelata",
                        "Mastophorus_muris", "Heterakis_spumosa", "Mesocestoides", 
                        "Catenotaenia_pusilla", "Hymenolepis", "Oxyurids", "Mix_Syphacia_Aspiculuris",
                        "Heligmosomoides_polygurus", "Latitude", "Longitude", "density")]

WormsDF2 <- na.omit(melt(WormsDF2, id = c("Longitude", "Latitude", "density")))

WormsDF2$prevalence <- 0
WormsDF2$prevalence[WormsDF2$value > 0] <- 1

aggworms2 <- aggregate(x = WormsDF2["prevalence"],
          by = WormsDF2[c("Latitude", "Longitude", "density", "variable")], FUN = function(x){sum(x)/length(x)})

ggplot(data=aggworms2, aes(x = density, y=prevalence, col = variable)) +
  geom_smooth(alpha = 0.1, size = 2) +
  geom_point(col = "black") +
  theme_classic() 

#***************************
## To finish

# link eimeria and prevalence 