library(httr)
library(RCurl)
library(dplyr)
library(ggplot2)
require(dplyr)
library(tidyverse)
library(reshape2)


#load in RT-qPCR data
RT <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/HZ18_RT-qPCR.csv"))
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

# load, process and add melting curve detection
MC <- read.csv(text = getURL("https://github.com/derele/Mouse_Eimeria_Databasing/blob/master/data/Eimeria_detection/HZ18_qPCR.csv"))

names(MC)[names(MC) == "Name"] <- "Mouse_ID"
MC <- MC %>% separate(Mouse_ID, c("CEWE", "AA", "Mouse_ID"))
MC$Mouse_ID <- sub("^", "AA_0", MC$Mouse_ID )
MC$CEWE <- NULL
MC$AA <- NULL
HZ18 <- merge(HZ18, MC, by = "Mouse_ID")
#subset MC 
Eim <- MC[,c(1,4,11)]
names(Eim)[names(Eim) == "Eimeria.presence.in.Caecum"] <- "inf"
RT.long <- merge(RT.long, Eim)
NE <- subset(x = RT.long$NE, subset = TRUE)
HZ18 <- merge(HZ18, RT.long)

# write out on Win home
write.csv(RT.long, file = "~/Mouse_Eimeria_Databasing/data/Gene_expression/HZ18_RT-qPCR_RTlong.csv", row.names = FALSE)
write.csv(HZ18, file = "~/Mouse_Eimeria_Databasing/data/Gene_expression/HZ18_complete.csv", row.names = FALSE)
#write out on Deb work
write.csv(RT.long, file = "~/Documents/Mouse_Eimeria_Databasing/data/Gene_expression/HZ18_RT-qPCR_RTlong.csv", row.names = FALSE)

ggplot(RT.long, aes(HI, NE, color = inf)) +
  geom_point() +
  geom_smooth() +
  facet_wrap("Target")

ggplot(RT.long, aes(deltaCtMmE_tissue, NE, color = inf))  +
  geom_point() + 
  facet_wrap("Target")+ 
  geom_smooth(method = "lm")

ggplot(RT.long, 
       aes(x = inf, y = NE, color = inf)) +
  geom_jitter() +
  geom_boxplot() +
  facet_wrap("Target", scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("HZ18_gene_downregulation")

############################################################### not used anymore
# eff.factor <- 1.9
# 
# RT.eff <-  eff.factor^(RT.wide[, c(refGenes, targetGenes)] * -1)
# 
# normIDX <- apply(RT.eff[, refGenes], 1, prod)^
#   (1/length(refGenes))
# 
# RT.norm <- RT.eff[, targetGenes] / normIDX
# 
# names(RT.norm) <- gsub("RT.Ct", "NE", names(RT.norm))
# 
# ## dropping everything but IDs and normalized values... look into SDs,
# ## non-normalized etc... if needed!!
# RT.norm <- cbind(Mouse_ID=RT.wide[, "Mouse_ID"], RT.norm)
# 
# # remove arbitrary RT columns, group, merge with HZ18
# HZ18 <- distinct(HZ18)
# HZ18 <- merge(HZ18, RT.norm, by = "Mouse_ID")
# #reshape RT.norm to graph
# RT.norm.long <- reshape(RT.norm, direction = "long", varying = c("NE.CXCR3", "NE.GBP2", "NE.IL-12b", "NE.IL-6", "NE.IRG6"), 
#                         idvar = "Mouse_ID")
# # RT.norm.long <- merge(RT.norm.long, detect, by = "Mouse_ID")
# RT.norm.long <- merge(RT.norm.long, HI, by = "Mouse_ID")
# names(RT.norm.long)[names(RT.norm.long)  == "time"] <- "Target"
# #graph (mus - eim = delta) (obsolete as we have melting curves)
# ggplot(data = RT.norm.long, aes(x = HI, y = NE, color = Target)) +
#   geom_point()
######################################################################################################


# Use MC in raw RT graphing
RT <- merge(RT, MC, by = 'Mouse_ID')
ggplot(data = RT, aes(x = HI, y = RT.Ct, color = Eimeria.presence.in.Caecum)) +
  geom_point() +
  geom_smooth() +
  labs(title = ("Raw CTs")) +
  facet_wrap("Target")

# graph Eim-mouse vs HI like Alice's bananas
MC <- merge(MC, HI, by = "Mouse_ID")
ggplot(MC, aes(HI, deltaCtMmE_tissue)) + 
  geom_point() +
  geom_smooth()
