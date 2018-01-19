## Load data from oocysts counting 
enas2015 <- read.csv("../raw_data/Eimeria_detection/Eimeria_oocysts_2015_Enas.csv")
lorenzo2015 <- read.csv("../raw_data/Eimeria_detection/Eimeria_oocysts_2015_Lorenzo.csv")
enas2016 <- read.csv("../raw_data/Eimeria_detection/Eimeria_oocysts_2016_part1_Enas.csv")
phuong2016 <- read.csv("../raw_data/Eimeria_detection/Eimeria_oocysts_2016_part2_Phuong.csv")
lorenzo2017 <- read.csv("../raw_data/Eimeria_detection/Eimeria_oocysts_2017_Lorenzo.csv")

options(scipen = 999)

eimeria_summary_df <- rbind(
  data.frame(Mouse_ID = enas2015$Mouse_ID,
             N_oocysts_in_feces = enas2015$oocysts_in_feces,
             OPG = enas2015$OPG,
             counter = enas2015$count, 
             year = enas2015$year),
  data.frame(Mouse_ID = lorenzo2015$Mouse_ID,
             N_oocysts_in_feces = lorenzo2015$oocysts_in_feces,
             OPG = lorenzo2015$OPG,
             counter = lorenzo2015$counter, 
             year = lorenzo2015$year),
  data.frame(Mouse_ID = enas2016$Mouse_ID, 
             N_oocysts_in_feces = enas2016$oocysts_in_feces,
             OPG = enas2016$OPG_if_0.4g, 
             counter =enas2016$count,
             year = enas2016$year),
  data.frame(Mouse_ID = phuong2016$Mouse_ID, 
             N_oocysts_in_feces = phuong2016$oocysts_in_feces,
             OPG = phuong2016$OPG_if_0.4g, 
             counter =phuong2016$count,
             year = phuong2016$year),
  data.frame(Mouse_ID = lorenzo2017$Mouse_ID, 
             N_oocysts_in_feces = lorenzo2017$oocysts_in_feces,
             OPG = lorenzo2017$OPG, 
             counter =lorenzo2017$counter,
             year = lorenzo2017$year))

# clean
eimeria_summary_df$year <- factor(eimeria_summary_df$year)
eimeria_summary_df <- eimeria_summary_df[!is.na(eimeria_summary_df$OPG),]

# write out
write.csv(x = eimeria_summary_df, 
          file = "../raw_data/Eimeria_detection/ALL_Eimeria_oocysts_2015_2016_2017.csv", 
          row.names = F)

# merge
a <- merge(data.frame(Mouse_ID = enas2015$Mouse_ID,
                      OPG.enas = enas2015$OPG,
                      year = enas2015$year),
           data.frame(Mouse_ID = lorenzo2015$Mouse_ID,
                      OPG.lorenzo = lorenzo2015$OPG,
                      year = lorenzo2015$year),
           by = c("Mouse_ID", "year"), all = TRUE)

a <- merge(a, data.frame(Mouse_ID = enas2016$Mouse_ID, 
                         OPG.enas = enas2016$OPG_if_0.4g, 
                         year = enas2016$year),
           by = c("Mouse_ID", "year", "OPG.enas"), all = TRUE)

a <- merge(a, data.frame(Mouse_ID = phuong2016$Mouse_ID, 
                         OPG.phuong = phuong2016$OPG_if_0.4g, 
                         year = phuong2016$year),
           by = c("Mouse_ID", "year"), all = TRUE)

a <- merge(a, data.frame(Mouse_ID = lorenzo2017$Mouse_ID, 
                         OPG.lorenzo = lorenzo2017$OPG,
                         year = lorenzo2017$year),
           by = c("Mouse_ID", "year", "OPG.lorenzo"), all = TRUE)

oocyst_summary_df <- a

## Load data from PCR
PCRdata <- read.csv("../raw_data/Inventory_contents_all.csv")

a <- PCRdata[c("X3_ID_mouse", "X13_18S_Seq", "X14_COI_Seq", "X16_ORF470_Seq")]
a[is.na(a)] <- 0

a$PCRpos <- "positive"
a$PCRpos[which(rowSums(a[2:4]) == 0)] <- "negative"

PCR_summary_df <- a
names(PCR_summary_df)[names(PCR_summary_df) == "X3_ID_mouse"] <- "Mouse_ID"

## merge all info (oocysts counting + PCR)
eimeria_detect <- merge(x = oocyst_summary_df, y = PCR_summary_df, 
                        by = "Mouse_ID", all = TRUE)

# calculate prevalence of different methods
prev <- function(x){table(x)[2]/sum(table(x))*100}

by(data = eimeria_detect$PCRpos, 
   INDICES = eimeria_detect$year,
   FUN = prev)

by(data = eimeria_detect$OPG.lorenzo > 0, 
   INDICES = eimeria_detect$year,
   FUN = prev)

by(data = eimeria_detect$OPG.enas > 0, 
   INDICES = eimeria_detect$year,
   FUN = prev)

by(data = eimeria_detect$OPG.phuong > 0, 
   INDICES = eimeria_detect$year,
   FUN = prev)

# venn diagram
# source("http://www.bioconductor.org/biocLite.R")
# biocLite("limma")
library(limma)

par(mfrow=c(2,2))

c1 <- cbind(lorenzo = eimeria_detect$OPG.lorenzo > 0,
            enas = eimeria_detect$OPG.enas > 0,
            PCR = eimeria_detect$PCRpos == "positive")

a <- vennCounts(c1)
a
vennDiagram(a)

## Lorenzo vs PCR
c2 <- cbind(lorenzo = eimeria_detect$OPG.lorenzo > 0,
            PCR = eimeria_detect$PCRpos == "positive")

a <- vennCounts(c2)
a
vennDiagram(a)

## Enas vs PCR
c3 <- cbind(enas = eimeria_detect$OPG.enas > 0,
            PCR = eimeria_detect$PCRpos == "positive")

a <- vennCounts(c3)
a
vennDiagram(a)
