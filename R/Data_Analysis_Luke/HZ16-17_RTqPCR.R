library(httr)
library(RCurl)
library(Rmisc)
library(tidyverse)
library(purrr)
library(ggplot2)
library(reshape2)

# load in runs
RT1 <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Gene_expression/HZ16-17_RT-qPCRs/HZ16-17_RT-qPCR1.CSV"
RT1 <- read.csv(text = getURL(RT1))

RT2 <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Gene_expression/HZ16-17_RT-qPCRs/HZ16-17_RT-qPCR2.CSV"
RT2 <- read.csv(text = getURL(RT2))

RT3 <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Gene_expression/HZ16-17_RT-qPCRs/HZ16-17_RT-qPCR3.CSV"
RT3 <- read.csv(text = getURL(RT3))

RT4 <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Gene_expression/HZ16-17_RT-qPCRs/HZ16-17_RT-qPCR4.CSV"
RT4 <- read.csv(text = getURL(RT4))

RT5 <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Gene_expression/HZ16-17_RT-qPCRs/HZ16-17_RT-qPCR5.CSV"
RT5 <- read.csv(text = getURL(RT5))

RT6 <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Gene_expression/HZ16-17_RT-qPCRs/HZ16-17_RT-qPCR6.CSV"
RT6 <- read.csv(text = getURL(RT6))

RT7 <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Gene_expression/HZ16-17_RT-qPCRs/HZ16-17_RT-qPCR7.CSV"
RT7 <- read.csv(text = getURL(RT7))

RT8 <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Gene_expression/HZ16-17_RT-qPCRs/HZ16-17_RT-qPCR8.CSV"
RT8 <- read.csv(text = getURL(RT8))

RT9 <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Gene_expression/HZ16-17_RT-qPCRs/HZ16-17_RT-qPCR9.CSV"
RT9 <- read.csv(text = getURL(RT9))
  
RT10 <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Gene_expression/HZ16-17_RT-qPCRs/HZ16-17_RT-qPCR10.CSV" 
RT10 <- read.csv(text = getURL(RT10))  
  
RT11 <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Gene_expression/HZ16-17_RT-qPCRs/HZ16-17_RT-qPCR11.CSV" 
RT11 <- read.csv(text = getURL(RT11)) 

# clean and merge it all
RT1$Ct.Mean.SYBR <- NULL
RT1$Ct.Dev..SYBR <- NULL
RT4$Ct.Mean.SYBR <- NULL
RT4$Ct.Dev..SYBR <- NULL
RT1$Amount.SYBR <- NULL
RT2$Amount.SYBR <- NULL
RT3$Amount.SYBR <- NULL
RT4$Amount.SYBR <- NULL
RT5$Amount.SYBR <- NULL
RT6$Amount.SYBR <- NULL
RT7$Amount.SYBR <- NULL
RT8$Amount.SYBR <- NULL
RT9$Amount.SYBR <- NULL
RT10$Amount.SYBR <- NULL
RT11$Amount.SYBR <- NULL

RT1$Pos <- NULL
RT2$Pos <- NULL
RT3$Pos <- NULL
RT4$Pos <- NULL
RT5$Pos <- NULL
RT6$Pos <- NULL
RT7$Pos <- NULL
RT8$Pos <- NULL
RT9$Pos <- NULL
RT10$Pos <- NULL
RT11$Pos <- NULL

RT <- rbind(RT1, RT2)
RT <- rbind(RT, RT3)
RT <- rbind(RT, RT4)
RT <- rbind(RT, RT5)
RT <- rbind(RT, RT6)
RT <- rbind(RT, RT7)
RT <- rbind(RT, RT8)
RT <- rbind(RT, RT9)
RT <- rbind(RT, RT10)
RT <- rbind(RT, RT11)

#remove negative controls
RT <- RT[!grepl("IRG6A", RT$Name),]
RT <- RT[!grepl("CXCR3", RT$Name),]
RT <- RT[!grepl("IL-12rb1", RT$Name),]
RT <- RT[!grepl("B-actin", RT$Name),]
RT <- RT[!grepl("GAPDH", RT$Name),]

# name columns to match other data sets (Mouse_ID, Target) + make RT.CT numeric
names(RT)[names(RT) == "Target.SYBR"] <- "Target"
names(RT)[names(RT) == "Ct.SYBR"] <- "RT.Ct"
names(RT)[names(RT) == "Name"] <- "Mouse_ID"
RT$RT.Ct <- as.numeric(as.character(RT$RT.Ct))

# remove NAs, calculate averages + save long
RT <- na.omit(RT)
RT.long <- RT %>% dplyr::group_by(Mouse_ID, Target) %>% dplyr::summarise(RT.Ct = mean(RT.Ct))
RT.long <- data.frame(RT.long)
RT.wide <- reshape(RT.long[, c("Target", "Mouse_ID","RT.Ct")],
                   timevar = "Target", idvar = "Mouse_ID", direction = "wide")
# set ref and target genes
refGenes <- c("RT.Ct.B-actin", "RT.Ct.GAPDH")
targetGenes <- c("RT.Ct.CXCR3", "RT.Ct.IL-12", "RT.Ct.IRG6")
# calculate ref genes in new column and subtract targets from HKG average, create new columns
require(dplyr)
RT.wide <- RT.wide %>% mutate(refMean = rowMeans(na.rm = TRUE, dplyr::select(RT.wide, refGenes)))
RT.wide <- data.frame(RT.wide)
refMean <- as.numeric(RT.wide$refMean)

# continue with averaging refgenes and subtracting targets from them
RT.wide$CXCR3 <- (RT.wide$refMean - RT.wide$RT.Ct.CXCR3)
RT.wide$IRG6 <- (RT.wide$refMean - RT.wide$RT.Ct.IRG6)
RT.wide$IL.12 <- (RT.wide$refMean - RT.wide$RT.Ct.IL.12)
# remove non normalized expressions
RT.wide$RT.Ct.CXCR3 <- NULL
RT.wide$RT.Ct.IRG6 <- NULL
RT.wide$RT.Ct.IL.12 <- NULL
RT.wide$RT.Ct.GAPDH <- NULL
RT.wide$RT.Ct.B.actin <- NULL
RT.wide$refMean <- NULL

RT.long <- melt(RT.wide, id.vars = "Mouse_ID")
names(RT.long)[names(RT.long) == "variable"] <- "Target"
names(RT.long)[names(RT.long) == "value"] <- "NE"

ggplot(RT.long, aes(x = NE, y = Mouse_ID)) +
  geom_point() +
  facet_wrap("Target")
##################################################################################

!!!!! start fix: fx by adding HZ18 genotypes to "HZgenotype" and dissection data from 2014-2017 to "diss"

#load in genotypes (make for 16-17)
genotypeURL <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Field_data/HZ10_HZ17_Genotypes_complete.csv"
HZgenotype <- read.csv(text = getURL(genotypeURL))
# subest by HI
HImus <- select(HZgenotype, HI, Mouse_ID)
#load in dissections
dissectionURL <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Field_data/HZ18_Dissections.csv"
HZ18dissection <- read.csv(text = getURL(dissectionURL))
#subset by columns relevant for mapping and gene expression
diss <- select(HZ18dissection, Mouse_ID, Latitude, Longitude, Sex, Status, Body_weight, Spleen, ASP, SYP, HET, MART, CP, HD, HM, MM, TM)
# merge HImus and diss
HImus <- merge(HImus, diss, by = "Mouse_ID")

!!!!!!! end fix

# correct names + add sample column
RT <- RT %>% separate(Name, c("CEWE", "AA", "Mouse_ID"))
RT$Mouse_ID <- sub("^", "AA_0", RT$Mouse_ID )
RT$AA <- NULL
names(RT)[names(RT) == "CEWE"] <- "tissue"
# calculate averages
RT <- RT %>% dplyr::group_by(Mouse_ID, Target.SYBR) %>% dplyr::summarise(Ct.SYBR = mean(Ct.SYBR))
#rename columns to merge by Mouse_ID
names(RT)[names(RT) == "Target.SYBR"] <- "Target"
names(RT)[names(RT) == "Ct.SYBR"] <- "RT.Ct"
# merge HImus and RT
HZ18 <- merge(RT, HImus)

# load in sample list


### Check against MC analysis (Lorenzo)
LorenzoMC <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/HZ16-17_InfInt_MC_Lorenzo%26Mert.csv"
LorenzoMC <- read.csv(text = getURL(LorenzoMC))

LorenzoMC <- merge(LorenzoMC, Positive)
TruePositives <- subset(LorenzoMC, Caecum == "pos")
write.csv(TruePositives, file = "~/Mouse_Eimeria_Databasing/Mouse_Eimeria_Databasing/data/Eimeria_detection/MC_verified_positives.csv")

############################## Add oocyst data
oocysts <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/Eimeria_oocysts_2015%262017_Lorenzo.csv"
oocysts <- read.csv(text = getURL(oocysts))
oocysts <- select(oocysts, Mouse_ID, mean_neubauer, OPG)
HZ16and17 <- merge(oocysts, HZ16and17, by = "Mouse_ID", all.y = TRUE)
# check for oocyst positive and qPCR negative flukes (keeps reducing the overall number... maybe stick to oocyst and HZ merge)
DoublePositive <- dplyr::filter(HZ16and17, qPCRstatus == "positive", OPG > 0)
DoublePositive$positive <- "double"
OocystPositive <- dplyr::filter(HZ16and17, qPCRstatus == "negative", OPG > 0)
OocystPositive$positive <- "oocyst"
qPCRPositive <- dplyr::filter(HZ16and17, qPCRstatus == "positive", OPG == 0)
qPCRPositive$positive <- "qPCR"
PositiveInvestigate <- rbind(DoublePositive, OocystPositive, qPCRPositive)
PositiveInvestigate <- distinct(PositiveInvestigate)
HZ16and17 <- merge(HZ16and17, PositiveInvestigate, by = "Mouse_ID")

################### Add melting cuve analysis data #######################

