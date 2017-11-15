## November 2017, Alice Balard
library(ggplot2)
library(ggmap)
library(ggrepel)

source("HMHZ_Functions.R")

MiceTable <- read.csv("../raw_data/MiceTable_2014to2017.csv")
TrapTable <- read.csv("../raw_data/TrapTable_2014to2017.csv")

## Let's split the localities VISITED and the localities where MICE WERE TRAPPED (before 2016)

## Cluster localities : how many times each LOCALITY was sampled?

# pairsloc <- pairwise.cluster.loc(trapsTOT) # errors with this model, to discuss
# loctime <- data.frame(table(apply(pairsloc, 2, sum)))
# loctime$Var1 <- as.numeric(as.character(loctime$Var1)) + 1 
# names(loctime) <- c("Repeatition over the years", "# localities")

measure(lon1 = 13.04390, lat1 = 51.93000, lon2 = 13.64000, lat2 = 52.71690) - measure(lon1 = 13.04, lat1 = 51.93, lon2 = 13.64, lat2 = 52.72)
measure(lon1 = 13.6761, lat1 = 52.4978, lon2 = 13.68, lat2 = 52.50)

# rounded to about 700 meters
TrapTable$Longitude_round_large <- round(TrapTable$Longitude, 2)
TrapTable$Latitude_round_large <- round(TrapTable$Latitude, 2)

# rounded to about 20 meters
TrapTable$Longitude_round_small <- round(TrapTable$Longitude, 3)
TrapTable$Latitude_round_small <- round(TrapTable$Latitude, 3)

nrow(unique(TrapTable[c("Latitude", "Longitude", "Year")]))
nrow(unique(TrapTable[c("Latitude_round_small", "Longitude_round_small", "Year")]))
nrow(unique(TrapTable[c("Latitude_round_large", "Longitude_round_large", "Year")]))

## We decide to round to 700 meters, and to consider the localities as in a cluster
countattemps <- TrapTable[c("Latitude_round_large", "Longitude_round_large", "Year")]
countattemps$attemps <- 1  

# how many time in one given year we sample twice at the same locality :
aggdata <- aggregate(x = countattemps$attemps, by = countattemps[c("Latitude_round_large", "Longitude_round_large", "Year")], FUN = sum)

# now we count only 1 trapping per year:
countattemps <- unique(countattemps)

# how many localities sampled per year? 
table(countattemps$Year)
data.frame(table(countattemps$Year))

ggplot(data.frame(table(countattemps$Year)), aes(x = Var1, y = Freq, fill = Freq), cex = 5) +
  geom_bar(stat = "identity", show.legend=F) +
  geom_text(aes(label = Freq), position = position_dodge(width=0.9), vjust=-0.25, size =10) +
  theme_classic(base_size = 18) +
  theme(axis.title.y=element_blank(), axis.title.x=element_blank()) +
  ggtitle(label = "Number of localities visited per year")











table(aggdata$attemps)

ggplot(data = data.frame(table(aggdata$Nmice)), aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "blue") +
  theme_black()

ggplot(data = data.frame(table(aggdata$Nmice)), aes(x = Var1, y = log10(Freq))) +
  geom_bar(stat = "identity", fill="blue" ) +
  theme_black()

table(aggdata$Nmice)




ggplot(loctime, aes(x = `Repeatition over the years`, y = `# localities`, fill = `# localities`)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label=Number), position=position_dodge(width=0.9), vjust=-0.25)
  theme_classic()


nrow(pairsloc)
ncol(pairsloc)

# 1. Sample size over the years: MICE and LOCALITIES

## Do the same for "tested for Eimeria" (cf Julia/Jenny/Phuong datasets)

## Later step: define the clusters:
## table(pairwise.cluster.loc(d = traps17))

## First step: with exact same GPS coordinates:
datacatch <- rbind(data.frame(table(paste(round(dissec14$X_Map,3), round(dissec14$Y_Map,3))), Year = 2014),
                   data.frame(table(paste(round(dissec15$X_Map,3), round(dissec15$Y_Map,3))), Year = 2015),
                   data.frame(table(paste(round(dissec16$Longitude,3), round(dissec16$Latitude,3))), Year = 2016),
                   data.frame(table(paste(round(dissec17$Longitude,3), round(dissec17$Latitude,3))), Year = 2017))
datacatch$Nloc <- 1

ggplot(datacatch, aes(x = Year, y = cumsum(Freq))) +
  geom_line(col = "lightgreen", size = 3) + 
  theme_black() + 
  scale_y_continuous(name = "Cumulative sum of mice caught") 

ggplot(datacatch, aes(x = Year, y = cumsum(Nloc))) +
  geom_line(col = "lightblue", size = 3) + 
  theme_black() + 
  scale_y_continuous(name = "Cumulative sum of localities sampled")
#***************************

# 2. Trapping attempts (from 2017 on):
sum(na.omit(as.numeric(as.character(traps17$Number_mus_caught))))
sum(na.omit(as.numeric(as.character(traps17$Number_traps_set))))

sum(na.omit(as.numeric(as.character(traps17$Number_mus_caught)))) /
  sum(na.omit(as.numeric(as.character(traps17$Number_traps_set)))) *100

# What else?

#***************************

# 3. Nbr mice per locality over the years?
aggdata <- aggregate(x = trapsTOT["Nmice"],
                     by = trapsTOT[c("Latitude", "Longitude")],
                     FUN = sum)
table(aggdata$Nmice)

ggplot(data = data.frame(table(aggdata$Nmice)), aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "blue") +
  theme_black()

ggplot(data = data.frame(table(aggdata$Nmice)), aes(x = Var1, y = log10(Freq))) +
  geom_bar(stat = "identity", fill="blue" ) +
  theme_black()

table(aggdata$Nmice)

#***************************

# 4. Nbr trapping success per locality over the years?

#***************************

# 4. Map





# Function to create barplots for different variables
mybarplot <- function(myfac, mytitle){
  ## set the levels in order we want
  dissec17 <- within(dissec17, 
                     myfac <- factor(myfac, 
                                     levels=names(sort(table(myfac), 
                                                       decreasing=TRUE))))
  ggplot(data = dissec17) +
    geom_bar(width = 0.5, aes(x=factor(1), fill=myfac)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(name = mytitle) +
    theme(legend.title=element_blank(), axis.text.x = element_blank())
}
mybarplot(dissec17ALL$Sex, "Sex")

# Map from traps
margin <- 0.3
area <- get_map(location =
                  c(min(dissect17MUS$Longitude - margin),
                    min(dissect17MUS$Latitude - margin),
                    max(dissect17MUS$Longitude + margin),
                    max(dissect17MUS$Latitude + margin)),
                source = "stamen", maptype="toner-lite",
                zoom = 8)

#plot the map :
ggmap(area) +
  geom_point(data = dissect17MUS, shape = 21, size = 4, fill = "green",
             aes(Longitude, Latitude)) +
  theme(legend.text=element_text(size=20)) +
  guides(fill=guide_legend(title=""))
  
# N average of captured mice (trapping success) in failed farms/successful farms
sum(na.omit(as.numeric(as.character(traps17$Number_rodents_caught))))
sum(na.omit(as.numeric(as.character(traps17$Number_mus_caught))))
sum(na.omit(as.numeric(as.character(traps17$Number_traps_set))))

# Distribution Nmice per localities
traps17$Number_mus_caught

# Plot of weight (color by species)


############### Crap below that line
## Alice Balard
## August 2017

# Source
InvCont <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/raw_data/Inventory_contents_all.csv")
GenandLoc <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/output_data/genDF_august2017.csv")

# Same names
names(InvCont)[c(1,3,4)] <- c("Year", "Code", "Mouse_ID")

# Correct errors :
InvCont$Transect <- gsub(pattern = " ", replacement = "", x = InvCont$Transect)
InvCont$Code <- gsub(pattern = " ", replacement = "", x = InvCont$Code)
InvCont$Mouse_ID <- gsub(pattern = " ", replacement = "", x = InvCont$Mouse_ID)
InvCont$Code <- as.character(InvCont$Code)
InvCont$Code <- gsub(pattern = "E_", replacement = "", x = InvCont$Code)

# Merge
FullDF <- merge(InvCont, GenandLoc, by = c("Mouse_ID", "Year"), all.x = TRUE)
FullDF$Code <- FullDF$Code.y
FullDF[which(is.na(FullDF$Code)),]$Code <- FullDF[which(is.na(FullDF$Code)),]$Code.x
FullDF <- FullDF[-c(which(names(FullDF) == "Code.x"), which(names(FullDF) == "Code.y"))]

# Any missing data?
FullDF[which(is.na(FullDF$HI)),]$Mouse_ID

# Write out
#write.csv(x = data.frame(Mouse_ID = FullDF[which(is.na(FullDF$HI)),]$Mouse_ID,
#           Code = FullDF[which(is.na(FullDF$HI)),]$Code,
#           Transect = FullDF[which(is.na(FullDF$HI)),]$Transect,
#           Year = FullDF[which(is.na(FullDF$HI)),]$Year),
#          file = "../output_data/information_missing_from_inventory_content.csv",
#          row.names = FALSE)

#########################################
# Prevalence analyses

## Tested = at least 1 marker tested :
FullDF$Tested <- FALSE
FullDF[which(!is.na(FullDF$X12_Ap5_PCR) | !is.na(FullDF$X13_COI_PCR) | !is.na(FullDF$X13_COI_PCR)),]$Tested <- TRUE

## Infected = Positive for at least 1 marker :
FullDF$INF <- FALSE
FullDF[which(FullDF$X12_Ap5_PCR == TRUE | FullDF$X12_18S_PCR == TRUE | FullDF$X13_COI_PCR == TRUE | FullDF$X15_ORF470_PCR == TRUE),]$INF <- TRUE

# Use DF:
useDF <- FullDF[which(FullDF$Tested %in% TRUE),]

##########
## For all year, for all transects :
myprevDF <- function(df, year, transect, threshold){
  data <- subset(df, df$Year == year & df$Transect == transect)    
  ## Prevalence table
  DF <- data.frame(
    Year = unique(data$Year),
    Transect = unique(data$Transect),
    N_Mmd = nrow(subset(data, data$HI < threshold)),
    percent_Mmd_inf = nrow(subset(data, data$HI < threshold & data$INF == TRUE)) / nrow(subset(data, data$HI < threshold)) *100,
    N_Hybrids = nrow(subset(data, data$HI >= threshold & data$HI <= 1- threshold)),
    percent_Hybrids_inf = nrow(subset(data, data$HI >= threshold & data$HI <= 1 - threshold & data$INF == TRUE)) / nrow(subset(data, data$HI >= threshold & data$HI <= 1- threshold))*100,
    N_Mmm = nrow(subset(data, data$HI > 1 - threshold)),
    percent_Mmm_inf = nrow(subset(data, data$HI > 1 - threshold & data$INF == TRUE)) / nrow(subset(data, data$HI > 1 - threshold))*100
  )
  DF
}

A <- myprevDF(df = FullDF, year = 2015, transect = "HZ_BR", threshold = 0.1)
B <- myprevDF(df = useDF, year = 2015, transect = "HZ_BAV", threshold = 0.1)
C <- myprevDF(df = useDF, year = 2016, transect = "HZ_BR", threshold = 0.1)

myprevDF2 <- function(df, year, threshold){
  data <- subset(df, df$Year == year)    
  ## Prevalence table
  DF <- data.frame(
    Year = unique(data$Year),
    N_Mmd = nrow(subset(data, data$HI < threshold)),
    percent_Mmd_inf = nrow(subset(data, data$HI < threshold & data$INF == TRUE)),
    N_Hybrids = nrow(subset(data, data$HI >= threshold & data$HI <= 1- threshold)),
    percent_Hybrids_inf = nrow(subset(data, data$HI >= threshold & data$HI <= 1 - threshold & data$INF == TRUE)),
    N_Mmm = nrow(subset(data, data$HI > 1 - threshold)),
    percent_Mmm_inf = nrow(subset(data, data$HI > 1 - threshold & data$INF == TRUE))
  )
  DF
}
myprevDF2(FullDF, 2016, 0.1)

A
B# Draw waffles :
library(waffle)
mywaffle <- function(letter, mytitle){
  vals <- c(as.numeric(letter[3] - letter[4]), as.numeric(letter[4]), 
            as.numeric(letter[5] - letter[6]), as.numeric(letter[6]),
            as.numeric(letter[7] - letter[8]), as.numeric(letter[8]))
  
  val_names <- sprintf("%s (%s)", c( "Mmd uninf", "Mmd inf", "Hybrid uninf", "Hybrid inf", "Mmm uninf", "Mmm inf"),
                       scales::percent(round(vals/sum(vals), 2)))
  names(vals) <- val_names
  waffle::waffle(vals, colors = c("dodgerblue4", "deepskyblue", "darkorchid4", "mediumorchid", "firebrick4", "firebrick1"),
                 title = mytitle)
}

mywaffle(A, "Brandenburg, 2015")
mywaffle(B, "Bavaria, 2015")
mywaffle(C, "Brandenburg, 2016")

summary(lm(formula = INF ~ Transect * Year * HI, data = useDF))

## Map of the data
source("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/R/HMHZ_Functions.R")

buildmap()
