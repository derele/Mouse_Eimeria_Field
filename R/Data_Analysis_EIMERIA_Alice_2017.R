## Load data from oocysts counting 
enas2015 <- read.csv("../raw_data/Eimeria_detection/Eimeria_oocysts_2015_Enas.csv")
lorenzo2015 <- read.csv("../raw_data/Eimeria_detection/Eimeria_oocysts_2015_Lorenzo.csv")
enas2016 <- read.csv("../raw_data/Eimeria_detection/Eimeria_oocysts_2016_part1_Enas.csv")
phuong2016 <- read.csv("../raw_data/Eimeria_detection/Eimeria_oocysts_2016_part2_Phuong.csv")
lorenzo2017 <- read.csv("../raw_data/Eimeria_detection/Eimeria_oocysts_2017_Lorenzo.csv")
alice <- read.csv("../raw_data/Eimeria_detection/Alice_newdilution_oocysts_counts_jan2018.csv")
alice <- na.omit(alice)
alice$alice$Year

options(scipen = 999)

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

a <- merge(a, data.frame(Mouse_ID = alice$Mouse_ID, 
                         OPG.alice = alice$OPG,
                         year = alice$year),
           by = c("Mouse_ID", "year"), all = TRUE)

oocyst_summary_df <- a
# correct name mistake!!
oocyst_summary_df$Mouse_ID <- gsub(" ", "", as.character(oocyst_summary_df$Mouse_ID))

## Load data from PCR
PCRdata <- read.csv("../raw_data/Eimeria_detection/Inventory_contents_all.csv")
# correct name mistake!!
PCRdata$X3_ID_mouse <- gsub(" ", "", as.character(PCRdata$X3_ID_mouse))

a <- PCRdata[c("X12_Ap5_PCR","X3_ID_mouse", "X13_18S_Seq", "X14_COI_Seq", "X16_ORF470_Seq")]

# Positive if 3 markers obtain sequences
a$PCRpos3markers <- "positive"
a$PCRpos3markers[which(
  rowSums(
    a[c("X13_18S_Seq", "X14_COI_Seq", "X16_ORF470_Seq")], 
    na.rm = T) == 0)] <- "negative"
a$PCRpos3markers[which(is.na(a$X13_18S_Seq) & 
                         is.na(a$X14_COI_Seq) & 
                         is.na(a$X16_ORF470_Seq))] <- NA


PCR_summary_df <- a[!is.na(a$PCRpos3markers),]
names(PCR_summary_df)[names(PCR_summary_df) == "X3_ID_mouse"] <- "Mouse_ID"

## merge all info (oocysts counting + PCR)
eimeria_detect <- merge(x = oocyst_summary_df, y = PCR_summary_df, 
                        by = "Mouse_ID", all = TRUE)

# venn diagram
# source("http://www.bioconductor.org/biocLite.R")
# biocLite("limma")
library(limma)

c1 <- cbind(lorenzo = eimeria_detect$OPG.lorenzo > 0,
            enas = eimeria_detect$OPG.enas > 0,
            alice = eimeria_detect$OPG.alice > 0,
            PCR = eimeria_detect$PCRpos3markers == "positive")

a <- vennCounts(c1)
a
vennDiagram(a, circle.col = 1:4, lwd = 3)

c2 <- cbind(lorenzo = eimeria_detect$OPG.lorenzo > 0,
            alice = eimeria_detect$OPG.alice > 0,
            PCR = eimeria_detect$PCRpos3markers == "positive")

a <- vennCounts(c2)
a
vennDiagram(a, circle.col = 1:4, lwd = 3)

############# Which samples to test? #############

# Oocysts counting found positive // PCR negative or not done:
eimeria2015 <- eimeria_detect[eimeria_detect$year == 2015,]

A <- eimeria2015[
  which(eimeria2015$PCRpos3markers == "negative" &
          (eimeria2015$OPG.alice != 0 | 
             eimeria2015$OPG.lorenzo != 0 |
             eimeria2015$OPG.enas != 0)),]

B <- eimeria2015[
  which(is.na(eimeria2015$PCRpos3markers) &
          (eimeria2015$OPG.alice != 0 | 
             eimeria2015$OPG.lorenzo != 0 |
             eimeria2015$OPG.enas != 0)),]

# Oocysts counting found negative // PCR positive or not done :
C <- eimeria2015[
  which(eimeria2015$PCRpos3markers == "positive" &
          (eimeria2015$OPG.alice == 0 | 
             eimeria2015$OPG.lorenzo == 0 |
             eimeria2015$OPG.enas == 0)),]

D <- eimeria2015[
  which(is.na(eimeria2015$PCRpos3markers) &
          (eimeria2015$OPG.alice == 0 | 
             eimeria2015$OPG.lorenzo == 0 |
             eimeria2015$OPG.enas == 0)),]

# No ocysts counting // PCR positive :
E <- eimeria2015[
  which(eimeria2015$PCRpos3markers == "positive" &
          (is.na(eimeria2015$OPG.alice == 0) | 
             is.na(eimeria2015$OPG.lorenzo == 0) |
             is.na(eimeria2015$OPG.enas == 0))),]

rbind(A, B,C, D,E)

write.csv(rbind(A, B,C, D,E), file = "../raw_data/Eimeria_detection/samples_to_test_temp.csv", row.names = F)

############# Comparison Lorenzo vs Alice counts #############
mydata <- na.omit(data.frame(L = eimeria_detect$OPG.lorenzo, 
                           A = eimeria_detect$OPG.alice,
                           mouse = eimeria_detect$Mouse_ID))
cor(mydata$L, mydata$A)

# http://www.sthda.com/english/wiki/correlation-test-between-two-variables-in-r

library("ggpubr")
ggscatter(mydata, x = "L", y = "A", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "lorenzo", ylab = "alice", size = 3)

library(reshape)
mydata <- melt(mydata, id = "mouse")

library(ggplot2)
ggplot(mydata, aes(x = mouse, y = value, col = variable)) +
  geom_point() +
  theme_classic()
