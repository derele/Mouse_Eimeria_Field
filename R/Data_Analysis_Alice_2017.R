## November 2017, Alice Balard
library(ggplot2)
library(ggmap)

traps17 <- read.csv(file = "../raw_data/HZ17_Mice_Trap.csv")
dissec17 <- read.csv(file = "../raw_data/HZ17_September_Mice_Dissection.csv")

# Check missing coordinates
data.frame(dissec17[is.na(dissec17$Latitude), ]$Capture, 
           dissec17[is.na(dissec17$Latitude), ]$Mouse_ID, 
           dissec17[is.na(dissec17$Latitude), ]$Address) 

# Number of mice, N other species
dissec17$host_type <- dissec17$Species
dissec17$host_type <- "other rodent"
dissec17 <- within(dissec17, host_type[Species == "Mus musculus"] <- "Mus musculus")

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

mybarplot(dissec17$host_type, "Host type")
table(dissec17$host_type)

## From here, work only with MUS MUSCULUS
dissect17MUS <- dissec17[dissec17$Species %in% "Mus musculus",]

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
  
# How many localities sampled?
length(unique(paste(traps17$Latitude, traps17$Longitude)))

# N average of captured mice (trapping success) in failed farms/successful farms

sum(na.omit(as.numeric(as.character(traps17$Number_rodents_caught))))
sum(na.omit(as.numeric(as.character(traps17$Number_mus_caught))))
sum(na.omit(as.numeric(as.character(traps17$Number_traps_set))))

# Distribution Nmice per localities
traps17$Number_mus_caught


# Plot of weight (color by species)




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






