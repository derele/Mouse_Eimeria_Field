library(httr)
library(RCurl)
library(dplyr)
library(ggplot2)
require(dplyr)
library(tidyverse)
library(reshape2)


#load in RT-qPCR data
RT <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Gene_expression/HZ18_RT-qPCR_RTlong.csv"))
# correct names + add sample column
RT <- RT %>% separate(Name, c("CEWE", "AA", "Mouse_ID"))
RT$Mouse_ID <- sub("^", "AA_0", RT$Mouse_ID )
RT$AA <- NULL
RT$CEWE <- NULL
# calculate averages
RT <- RT %>% dplyr::group_by(Mouse_ID, Target.SYBR) %>% dplyr::summarise(Ct.SYBR = mean(Ct.SYBR))
#rename columns to merge by Mouse_ID
names(RT)[names(RT) == "Target.SYBR"] <- "Target"
names(RT)[names(RT) == "Ct.SYBR"] <- "RT.Ct"


RT <- data.frame(RT)
# RT <- RT %>% drop_na(RT.Ct)
RT.wide <- reshape(RT[, c("Target", "Mouse_ID","RT.Ct")],
                   timevar = "Target", idvar = "Mouse_ID", direction = "wide")


# subtract ref genes and make long for graphing
refGenes <- c("RT.Ct.beta-Actin", "RT.Ct.GAPDH")
targetGenes <- c("RT.Ct.CXCR3", "RT.Ct.GBP2", "RT.Ct.IL.12b", 
                 "RT.Ct.IL.6", "RT.Ct.IRG6")
RT.wide <- RT.wide %>% mutate(refMean = rowMeans(select(., refGenes)))
RT.wide <- data.frame(RT.wide)
refMean <- as.numeric(RT.wide$refMean)

RT.wide$CXCR3 <- (RT.wide$refMean - RT.wide$RT.Ct.CXCR3)
RT.wide$IRG6 <- (RT.wide$refMean - RT.wide$RT.Ct.IRG6)
RT.wide$IL.12b <- (RT.wide$refMean - RT.wide$RT.Ct.IL.12b)
RT.wide$GBP2 <- (RT.wide$refMean - RT.wide$RT.Ct.GBP2)
RT.wide$IL.6 <- (RT.wide$refMean - RT.wide$RT.Ct.IL.6)
RT.wide$RT.Ct.beta.Actin <- NULL
RT.wide$RT.Ct.GAPDH <- NULL
RT.wide$refMean <- NULL
RT.wide$RT.Ct.CXCR3 <- NULL
RT.wide$RT.Ct.GBP2 <- NULL
RT.wide$RT.Ct.IL.12b <- NULL
RT.wide$RT.Ct.IL.6 <- NULL
RT.wide$RT.Ct.IRG6 <- NULL

names(RT.wide)[names(RT.wide) == "CXCR3"] <- "CXCR3"
names(RT.wide)[names(RT.wide) == "GBP2"] <- "GBP2"
names(RT.wide)[names(RT.wide) == "IL.12b"] <- "IL-12b"
names(RT.wide)[names(RT.wide) == "IL.6"] <- "IL-6"
names(RT.wide)[names(RT.wide) == "IRG6"] <- "IRG6"

RT.long <- melt(RT.wide, id.vars = "Mouse_ID")
names(RT.long)[names(RT.long) == "variable"] <- "Target"
names(RT.long)[names(RT.long) == "value"] <- "NE"

write.csv(RT.wide, "/Users/Luke Bednar/Mouse_Eimeria_Databasing/data/Gene_expression/HZ18_RT-qPCR_complete.csv")




ggscatter(RT, x = "deltaCtMmE_tissue", y = "NE", add = "reg.line") +
  facet_grid(inf~Target, scales = "free")+
  stat_cor(label.x =-5, label.y = 0) +
  stat_regline_equation(label.x = -5, label.y = 0) + 
  labs(y = "NE", x = "deltaCT = Mouse - Eimeria") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("expression vs delta")

ggscatter(RT, x = "inf", y = "NE", add = "reg.line") +
  facet_grid(~Target, scales = "free")+
  geom_boxplot() +
  stat_cor(label.x =-5, label.y = 0) +
  stat_regline_equation(label.x = -5, label.y = 0) + 
  labs(y = "NE", x = "deltaCT = Mouse - Eimeria") +
  theme(axis.text=element_text(size=12, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("Genes during infections")
