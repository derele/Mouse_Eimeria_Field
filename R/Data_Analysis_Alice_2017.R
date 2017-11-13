## November 2017, Alice Balard
source("HMHZ_Functions.R")

library(ggplot2)
library(ggmap)
library(ggrepel)

#*************************** Input data CLEAN ONCE FOR ALL
jardaHI14_16 <- read.csv("../raw_data/HIforEH_May2017.csv")
cleanloc <- read.csv("../output_data/all_clean_localities.csv") # delete and summarize



Dissection14to17 <- data.frame(Mouse_ID = NA, Transect = NA, Latitude = NA, Longitude = NA,
                               Sex = NA, Status = NA, Species = NA, Year = NA,
                               Capture_date = NA, Dissection_date = NA, Ectoparasites = NA,
Body_weight = NA, Body_length = NA, Tail_length = NA,

                               [19] "Spleen_mass"                        "Left_Testis_mass"                   "Right_Testis_mass"                 
                               [22] "Left.epididymis.weight"             "Seminal.vesicle.weight"             "Left.ovarium.weight"               
                               [25] "Right.ovarium.weight"               "Total"                              "Embryo_left"                       
                               [28] "Embryo_right"                       "Feces_weight"                       "Worms_presence"                    
                               [31] "Unknown_cecum"                      "Unknown_colon"                      "Unknown_SI"                        
                               [34] "Syphacia_cecum_colon"               "Cysticercus_tenia_teniformis_liver" "Catotenia_pusilla_SI"              
                               [37] "Mastophorus_stomach_SI"             "Heterakis_spumosa_colon_cecum"      "Tapeworm_SI"                       
                               [40] "Trichuris_cecum"                    "Aspiculuris_cecum"                  "Aspiculuris_colon"                 
                               [43] "Hymenolepis_SI"                     "Rodentolepis_liver_digtract"        "Mesocoides_body_cavities"          
                               [46] "Mesocoides_lungs"                   "Heigmosomoides_polyguis"            "H_diminita."                       
                               [49] "Notes"  
)

names(dissec17)

## Dissections
dissec14 <- read.csv("../raw_data/HZ14_Mice 31-12-14_dissections.csv")

  
  
  dissec14$Transect, dissec14$State))

Dissection14to17
## NB: Julia has 87 from HZ_BR (all here) + 165 from HZ_CZ 10 days (here 191 in 10 days...)
dissec15 <- read.csv("../raw_data/HZ15_Mice_Parasite.csv")
# gen15 <- read.csv("../raw_data/Genotypes_Bav2015.csv")
## NB: Jenny has 119 from HZ_BAV (here 276 +++ breeding) + 124 from HZ_BRA (here 152 some dissected in Studenec)
dissec16 <- read.csv("../raw_data/HZ16_Mice_18-07-16_dissections.csv")
dissec16 <- merge(dissec16, cleanloc, by = "Code", all.x = TRUE)
tocomplete <- jardaHI14_16[jardaHI14_16$PIN %in% dissec16[is.na(dissec16$Latitude),"ID_mouse"],c("PIN","Xmap", "Ymap")]
dissec16[dissec16$ID_mouse %in% tocomplete$PIN,"Latitude"] <- tocomplete$Ymap
dissec16[dissec16$ID_mouse %in% tocomplete$PIN,"Longitude"] <- tocomplete$Xmap
rm(tocomplete)

dissec17ALL <- read.csv(file = "../raw_data/HZ17_September_Mice_Dissection.csv")
dissec17 <- dissec17ALL[dissec17ALL$Species %in% "Mus musculus",]

names(dissec17)
names(dissec16)
names(dissec15)
names(dissec14)



## Traps
traps14to16 <- read.csv("../output_data/HZ14-16_localities.csv")[-1] ## REALLY WRONG UNIQUE HI
traps17 <- read.csv(file = "../raw_data/HZ17_Mice_Trap.csv")
trapsTOT <- data.frame(Year = c(rep(2017, nrow(traps17)), traps14to16$Year),
                       Latitude = c(traps17$Latitude, traps14to16$Latitude),
                       Longitude = c(traps17$Longitude, traps14to16$Longitude),
                       HI = c(rep(NA, nrow(traps17)), traps14to16$HI),
                       Nmice = c(traps17$Number_mus_caught, traps14to16$n.mice))
#***************************

# Cluster localities
pairsloc <- pairwise.cluster.loc(trapsTOT)

loctime <- data.frame(table(apply(pairsloc, 2, sum)))
loctime$Var1 <- as.numeric(as.character(loctime$Var1)) + 1 
names(loctime) <- c("Repeatition over the years", "# localities")

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
