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



 
