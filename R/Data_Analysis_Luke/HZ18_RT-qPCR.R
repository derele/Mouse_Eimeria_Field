library(httr)
library(RCurl)
library(dplyr)
library(ggplot2)
require(dplyr)

#load in genotypes
genotypeURL <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Field_data/HZ18_Genotypes.csv"
HZ18genotype <- read.csv(text = getURL(genotypeURL))
# subest by HI
HImus <- select(HZ18genotype, HI, Mouse_ID)
#load in dissections
dissectionURL <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Field_data/HZ18_Dissections.csv"
HZ18dissection <- read.csv(text = getURL(dissectionURL))
#subset by columns relevant for mapping and gene expression
diss <- select(HZ18dissection, Mouse_ID, Latitude, Longitude, Sex, Status, Body_weight, Spleen, ASP, SYP, HET, MART, CP, HD, HM, MM, TM)
# merge HImus and diss
HImus <- merge(HImus, diss, by = "Mouse_ID")
#load in RT-qPCR data
RTURL<- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/HZ18_RT-qPCR.csv"
RT <- read.csv(text = getURL(RTURL))
# correct names + add sample column
RT$Name <- sapply(RT$Name, as.character)
RT$Name <- gsub("(\\D{4})_\\D{2}_\\d{3}", "\\4", RT$Sample)
cells$tissue <- gsub("\\d+: (mLN|spleen)_(\\d{3})_\\d{3}.fcs", "\\1", cells$Sample)
