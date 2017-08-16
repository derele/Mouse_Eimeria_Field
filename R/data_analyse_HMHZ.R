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



Myprev <- function(data, HI){
  ## Positive for at least 1 marker :
  InvPos <- subset(data, X12_Ap5_PCR == TRUE | X12_18S_PCR == TRUE | X13_COI_PCR == TRUE | X15_ORF470_PCR == TRUE)
  ## Prevalence :
  
}

Myprev(subset(InvCont, X1_Year == 2015 & InvCont$Transect == "HZ_BR"))

Mytab <- function(myyear){
  data <- subset(InvCont, X1_Year == myyear)
  # Info
  N_mice_tot <- length(unique(data$X3_ID_mouse))
  N_localities <- length(unique(data$X2_Code))
  # For each category :
  mysubtab <- function(HI){
    data2 <- subset(data, )
    
  }
  
  
  ## Positive for at least 1 marker :
  InvPos <- subset(data, X12_Ap5_PCR == TRUE | X12_18S_PCR == TRUE | X13_COI_PCR == TRUE | X15_ORF470_PCR == TRUE)
  
  
  Year = data$X1_Year
  N_Mmd = 
    N_Mmd_inf
  N_Hybrids
  N_Hybrids_inf
  N_Mmm
  N_Mmm_inf
  
  