library(httr)
library(RCurl)
library(dplyr)
library(ggplot2)
require(dplyr)
library(tidyverse)
library(reshape2)

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
RT <- RT %>% dplyr::group_by(Mouse_ID, Target.SYBR) %>% dplyr::summarise(Ct.SYBR = mean(Ct.SYBR))
#rename columns to merge by Mouse_ID
names(RT)[names(RT) == "Target.SYBR"] <- "Target"
names(RT)[names(RT) == "Ct.SYBR"] <- "RT.Ct"
# merge HImus and RT
HZ18 <- merge(RT, HImus)
#start graphing
ggplot(data = HZ18, aes(x = HI, y = RT.Ct)) +
  geom_point() + 
  facet_wrap("Target")
# # add infection intensity data (Eimeria - Mouse)
# detectURL <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/Svenja/joined_qPCR_tables_cecum.csv"
# detect <- read.csv(text = getURL(detectURL))
# #cleanup and merge
# detect <- detect %>% separate(Mouse_ID, c("CEWE", "AA", "Mouse_ID"))
# detect$Mouse_ID <- sub("^", "AA_0", detect$Mouse_ID )
# names(detect)[names(detect) == "CEWE"] <- "tissue"
# detect$AA <- NULL
# detect$Ct.SYBR <- NULL
# detect[,4:8] <- NULL
# names(detect)[names(detect) == "Ct.Mean.SYBR"] <- "inf.Ct"
# names(detect)[names(detect) == "Ct.Dev..SYBR"] <- "inf.Ct.Dev"
# detect <- detect %>% drop_na(delta)
# HZ18 <- merge(detect, HZ18, by = "Mouse_ID")
# # graph
# ggplot(data = HZ18, aes(x = HI, y = RT.Ct, color = delta)) +
#   geom_point() + 
#   facet_wrap("Target")
# calculate endogenous controls and log expression
RT <- data.frame(RT)
RT <- RT %>% drop_na(RT.Ct)
RT.wide <- reshape(RT[, c("Target", "Mouse_ID","RT.Ct")],
                   timevar = "Target", idvar = "Mouse_ID", direction = "wide")
# side trip to look at raw CT values before transformation
HI <- HImus[, 1:2]
RT <- merge(RT, HI, by = 'Mouse_ID')
ggplot(RT, aes(HI, RT.Ct)) +
  geom_point() +
  geom_smooth() +
  facet_wrap('Target')
# subtract ref genes and make long for graphing
refGenes <- c("RT.Ct.beta-Actin", "RT.Ct.GAPDH")
targetGenes <- c("RT.Ct.CXCR3", "RT.Ct.GBP2", "RT.Ct.IL.12b", 
                 "RT.Ct.IL.6", "RT.Ct.IRG6")
RT.wide <- RT.wide %>% mutate(refMean = rowMeans(select(., refGenes)))
RT.wide <- data.frame(RT.wide)
refMean <- as.numeric(RT.wide$refMean)
RT.wide$RT.Ct.beta.Actin <- NULL
RT.wide$RT.Ct.GAPDH <- NULL
RT.wide[,2:6] <- sweep(RT.wide[2:6],1,refMean,'-')
RT.wide$refMean <- NULL

names(RT.wide)[names(RT.wide) == "RT.Ct.CXCR3"] <- "CXCR3"
names(RT.wide)[names(RT.wide) == "RT.Ct.GBP2"] <- "GBP2"
names(RT.wide)[names(RT.wide) == "RT.Ct.IL.12b"] <- "IL-12b"
names(RT.wide)[names(RT.wide) == "RT.Ct.IL.6"] <- "IL-6"
names(RT.wide)[names(RT.wide) == "RT.Ct.IRG6"] <- "IRG6"

RT.long <- melt(RT.wide, id.vars = "Mouse_ID")
names(RT.long)[names(RT.long) == "variable"] <- "Target"
names(RT.long)[names(RT.long) == "value"] <- "NE"
RT.long <- merge(RT.long, HI)
# load, process and add melting curve detection
MCURL <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Eimeria_detection/Svenja/table_ct_and_more.csv"
MC <- read.csv(text = getURL(MCURL))
names(MC)[names(MC) == "Name"] <- "Mouse_ID"
MC <- MC %>% separate(Mouse_ID, c("CEWE", "AA", "Mouse_ID"))
MC$Mouse_ID <- sub("^", "AA_0", MC$Mouse_ID )
MC$CEWE <- NULL
MC$AA <- NULL
HZ18 <- merge(HZ18, MC, by = "Mouse_ID")
#subset MC 
Eim <- MC[,c(1,4)]
names(Eim)[names(Eim) == "Eimeria.presence.in.Caecum"] <- "inf"
RT.long <- merge(RT.long, Eim)

# write out on Win home
write.csv(RT.long, file = "~/Mouse_Eimeria_Databasing/data/Gene_expression/HZ18_RT-qPCR_RTlong.csv", row.names = FALSE)

ggplot(RT.long, aes(HI, NE, color = inf)) +
  geom_point() +
  geom_smooth() +
  facet_wrap("Target")

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
