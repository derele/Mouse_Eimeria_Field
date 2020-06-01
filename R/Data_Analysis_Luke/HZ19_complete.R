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
FACS.long <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Field_data/HZ19_FACS_long_mln.csv"))
FACS$X <- NULL
FACS.long$X <- NULL

immuno <- merge(qPCR, ELISA_CEWE, all = T)
immuno <- merge(immuno, FACS, all = T)

write.csv(immuno, "/Users/Luke Bednar/Mouse_Eimeria_Databasing/data/HZ19_immuno.csv")

# turn tables into long and merge to graph

immuno.long <- merge(qPCR, ELISA_CEWE)
immuno.long <- merge(immuno.long, FACS.long)
write.csv(immuno.long, "/Users/Luke Bednar/Mouse_Eimeria_Databasing/data/wild_immuno_long.csv")
# let's have a look
ggscatter(immuno.long, x = "IFNy", y = "delta", add = "reg.line", color = "MC.Eimeria") +
  facet_wrap(~MC.Eimeria)+
  stat_cor(label.x = 50, label.y = 0) +
  stat_regline_equation(label.x = 50, label.y = 2) + 
  ggtitle("HZ19 infections vs IFNy")

# now FACS look at the delta vs populations
ggscatter(immuno.long, x = "delta", y = "counts", add = "reg.line", color = "MC.Eimeria") +
  facet_grid(~pop)+
  stat_cor(label.x = -20, label.y = 0) +
  stat_regline_equation(label.x = -20, label.y = 5) + 
  ggtitle("HZ19 infections vs IFNy")

ggscatter(subset(immuno.long, immuno.long$MC.Eimeria == "TRUE"), x = "counts", y = "delta", add = "reg.line", color = "MC.Eimeria") +
  facet_wrap(~pop)+
  stat_cor(label.x = -20, label.y = 0) +
  stat_regline_equation(label.x = -20, label.y = 5) + 
  ggtitle("HZ19 infections vs IFNy")

ggscatter(subset(immuno.long, immuno.long$MC.Eimeria == "FALSE"), x = "counts", y = "delta", add = "reg.line", color = "MC.Eimeria") +
  facet_wrap(~pop)+
  stat_cor(label.x = -20, label.y = 0) +
  stat_regline_equation(label.x = -20, label.y = 5) + 
  ggtitle("HZ19 infections vs IFNy")
