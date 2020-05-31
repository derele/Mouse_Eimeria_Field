library(Rmisc)
library(httr)
library(RCurl)
library(plyr)
library(dplyr)
library(reshape2)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(stats)
library(ggsignif)
library(ggpmisc)

qPCR <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/HZ19_qPCR.csv"))
qPCR$X <- NULL
# doesn"t exist yet
#RT <- read.csv(text = getURL(""))

ELISA_CEWE <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/ELISAs/HZ19_CEWE_ELISAs_complete.csv"))
ELISA_CEWE$X <- NULL

FACS <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Field_data/HZ19_FACS_complete.csv"))
FACS$X <- NULL

immuno <- merge(qPCR, ELISA_CEWE, all = T)
immuno <- merge(immuno, FACS, all = T)

write.csv(immuno, "/Users/Luke Bednar/Mouse_Eimeria_Databasing/data/HZ19_immuno.csv")

# turn tables into long and merge to graph

immuno.long <- merge(qPCR, ELISA_CEWE)
# let's have a look
ggscatter(immuno.long, x = "IFNy", y = "delta", add = "reg.line", color = "MC.Eimeria") +
  facet_wrap(~MC.Eimeria)+
  stat_cor(label.x = 50, label.y = 0) +
  stat_regline_equation(label.x = 50, label.y = 2) + 
  ggtitle("HZ19 infections vs IFNy")

# now let's add FACS and look at the delta vs populations
FACS.long <- melt(FACS,
               direction = "long",
               varying = list(names(FACS)[2:13]),
               v.names = "cell.pop",
               na.rm = T, value.name = "counts", 
               id.vars = c("Mouse_ID", "Position", "infHistory"))
