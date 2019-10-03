library(httr)
library(RCurl)
library(dplyr)
library(ggplot2)
require(dplyr)
library(tidyverse)

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
RT <- RT %>% separate(Name, c("CEWE", "AA", "Mouse_ID"))
RT$Mouse_ID <- sub("^", "AA_0", RT$Mouse_ID )
RT$AA <- NULL
names(RT)[names(RT) == "CEWE"] <- "tissue"
# calculate averages
RT <- data.frame(RT %>% group_by(Target.SYBR, Mouse_ID, add = TRUE) %>% 
                       summarize(SD = sd(Ct.SYBR),
                                 Ct.Mean = mean(Ct.SYBR)))
#rename columns to merge by Mouse_ID
names(RT)[names(RT) == "Target.SYBR"] <- "Target"
# merge HImus and RT
HZ18 <- merge(RT, HImus)
#start graphing
ggplot(data = HZ18, aes(x = HI, y = Ct.Mean)) +
  geom_point() + 
  facet_wrap("Target")
# add infection intensity data
