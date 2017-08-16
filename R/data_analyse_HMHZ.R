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

##########
## For all year, for all transects :
myprevDF <- function(year, transect){
  data <- subset(FullDF, FullDF$Year == year & FullDF$Transect == transect)    
  ## Prevalence table
  DF <- data.frame(
    Year = unique(data$Year),
    Transect = unique(data$Transect),
    N_Mmd = nrow(subset(data, data$HI < 0.2)),
    N_Mmd_inf = nrow(subset(data, data$HI < 0.2 & data$INF == TRUE)),
    N_Hybrids = nrow(subset(data, data$HI >= 0.2 & data$HI <= 0.8)),
    N_Hybrids_inf = nrow(subset(data, data$HI >= 0.2 & data$HI <= 0.8 & data$INF == TRUE)),
    N_Mmm = nrow(subset(data, data$HI > 0.8)),
    N_Mmm_inf = nrow(subset(data, data$HI > 0.8 & data$INF == TRUE))
  )
  DF
}

myprevDF(year = 2014, transect = "HZ_BR")
myprevDF(year = 2015, transect = "HZ_BR")
myprevDF(year = 2015, transect = "HZ_BAV")
myprevDF(year = 2016, transect = "HZ_BR")