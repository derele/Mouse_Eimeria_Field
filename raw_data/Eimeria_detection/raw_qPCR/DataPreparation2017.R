# I am just looking at the Tm values to see
# 
# Standart deviation calculated on the triplicates
# e.g. for 3 mouse triplicates on 1 plate it gives 3 Ct, 1 Ct mean and 1 Ct sd
# 
# emanuel [2:34 PM]
# Okay! Got it! But you don't need Eimeria data. No Eimeria data is good (unifected) data. Do yoiu set this at an arbitrary high ct?
# Negative Eimeria detection can also have a high *SD of the* ct. E.g.NA, 35, 40 -> perfectly negative sample. (edited)
# Bad melting <76 for the 35 sample -> even better, more evidence this is negative.
# 
# alice [2:36 PM]
# hum. but even bad you would expect the 3 triplicates of eimeria on 1 plate to have close Ct right? as it's technical replicates?
#   
#   emanuel [2:37 PM]
# No, because it's noise. Noise can have a high SD of the ct.
# 
# alice [2:37 PM]
# huuuuuuOOOOOOOOOOOOOOOOOOuuuuuuuu
# but that makes the choice of Ct_eimeria value horrible no?
# you just take the mean, whatever the sd, and accept that there is intense noise?
# 
# emanuel [2:38 PM]
# I'd predict that if you have something like 35, 39, no measurement. the 35 rep would have a melt ct peak of 75.
# 
# alice [2:38 PM]
# yes I saw that punctually
# (when looking at the pdf reports manually)
# 
# emanuel [2:39 PM]
# Yes on the above. Or you ~discard~ *set to ct=NA* melt ct <76 first then set NA to max (ct) or 40 or whatever. Then you take the mean. (edited)
# At the end all this data will be considerd noise and cancled out. qPCR can only be used to analyse positive.
# 
# alice [2:41 PM]
# okidoki. I will do that then!
#   
#   emanuel [2:42 PM]
# The discard melt ct <76 is harsh and questionable. Maybe you want to play with that. E.g. melt ct <76 and detection ct > 30 set to 40.
# "discard" above should likely be "set to 40".
# 
# alice [2:43 PM]
# one question: sometimes, one sample did not work on one plate, so the full plate was replicated, even samples that worked first time. So at the end we have up to 4 replications for some samples. Should I keep all, or just one per sample? I was going for the last, but then how should I select the replicate to peak? I thought first based on sd, but this dicsussion made me doubt. Could be either on lower mouse Ct sd, or on date (e.g. the first one that worked), or...
# melt ct <76 & detection ct > 30 => Ct set to 40
# 
# emanuel [2:44 PM]
# ~Yes.~ (edited)
# 
# alice [2:44 PM]
# melt ct >76 & detection ct < 30 => Ct set to 40 ?? (edited)
# 
# emanuel [2:45 PM]
# No. melt ct >76 is an indicator of proper amplification.
# 
# alice [2:45 PM]
# ok I just got lost there. But yes of course
# 
# emanuel [2:45 PM]
# detect ct of 32 could still be correct then.
# 
# alice [2:45 PM]
# I meant more : melt ct = NA & detection ct > 30 => Ct set to 40
# 
# emanuel [2:46 PM]
# ~Yes.~ No proper drop in the melting curve indicates likely no product. But you can also leave this ct as it is! (edited)
# 
# alice [2:46 PM]
# melt ct = NA & detection ct < 30 => :upside_down_face:  ?
#   
#   emanuel [2:48 PM]
# Exactly head spinning. Leave cts > 30 as they are, they are fine. They are the usual noise.
# I think I was wrong above. Why should one do that. That indicates more evidence for non-infection that we have.
# Sorry for changing my mind that quickly.
# Ahhh ... now I remember. My reasoning was based on combining with other technical replicates.
# Whatever you do, you have to do it by sample.
# or for the non-by crowd `for index in all_samples`
# so... for Eimeria amplification of each sample: (edited)
# a) do you have a tech replicate with detection ct? (edited)
# No -> set 40 be happy
# Yes -> use only those with ct and go to b) (edited)
# b) is the melting ct of any sample >76
# No -> if detection ct > 30 use them if detection ct <30 flag as problematic (are their any? I know this is not aproper flow).
# YES -> discard all others use those and go to c)
# 
# alice [3:00 PM]
# (by sample you mean triplicate? cause I often have within 1 triplicate 1 with Ct and not the others)
# (I would consider that not infected)
# 
# emanuel [3:00 PM]
# c) calculate mean and sd
# By mouse and tissue. By biological entity about which you want to figure out something.
# 
# alice [3:02 PM]
# (yes but in 1 plate, 1 mouse, 1 tissue, you have 3 values. Sometimes not identical/with or without melting curve e.g)
# 
# emanuel [3:03 PM]
# This stepwise progression puts more emphasis on any possible doubt about negativity. If you'd include 40s (or whatever) in a mean you would do this put more emphasis on indication of absence. (edited)
# Do you want to consider plate?
# Of course you have to discard "simply wrong" amplifications. E.g. those starting with too high template and leading to a dodgy ct detection curve...
# 
# alice [3:05 PM]
# ok, thanks for all that, I will do that now and see what we end up with...
# I let you know in few hours :wink:
# 
# emanuel [3:07 PM]
# Great! Just one thing: You wont be able to deal with straight visible error. Too high template concentration leads e.g. to slowly increasing fluorescece which might reach ct at very early cycle. It's wrong anyways. Should be discarded.
# That's another problem then the Eimeria background. You might address this second problem successfully with data analysis (instead of the arbitrary selections others make).
# 
# alice [3:09 PM]
# ok
# Just looked at some raw qPCR report with Victor. The melt CT >76 criterion for eimeria amplifications is clearly wrong.
# Most proper amplifications for positive (other assays) samples had melt CT ~75 and low detection CT. Nicely positive samples.
# I was wrong :face_with_monocle:
#   Something in your data prep might be wrong.
# Maybe we should rethink and discuss this together ... 


# import data
rawData <- read.csv("./LorenzoRAW/CSVFiles/TotalLorenzo.csv", stringsAsFactors = F,
                 na.strings = c("NA", "", " ", "-"))

##### Prepare data #####
# Column file name
names(rawData)[names(rawData) %in% "QPCR01.06.2018.XLS.csv"] <- "fileName"

# Annoying spaces
rawData[,names(rawData)  %in% "fileName"] <- gsub(" ", "", rawData[,names(rawData) %in% "fileName"])
rawData[,names(rawData)  %in% "Target.SYBR"] <- gsub(" ", "", rawData[,names(rawData) %in% "Target.SYBR"])

# Remove inside headers
rawData <- rawData[rawData$Name != "Name",]

rawData[rawData$fileName %in% "QPCR03.05.2018_complete.csv", "fileName"] <- "QPCR03.05.2018.XLS.csv"  

rawData$fileName[rawData$fileName %in% "QPCR311-317.XLS.csv"] <- "QPCR09.05.2018.XLS.csv"

## Add melting curves infos when available
rawMeltData <- read.csv("./LorenzoRAW/MeltingCurves/MeltingCurvesLorenzo.csv",
                        na.strings = c("NA", "", " ", "-"))

names(rawMeltData)[names(rawMeltData) %in% "QPCR01.06.2018_MeltCurve.csv"] <- "fileName"

# Annoying spaces
rawMeltData[,names(rawMeltData)  %in% "fileName"] <- gsub(" ", "", rawMeltData[,names(rawMeltData) %in% "fileName"])

# Remove inside headers
rawMeltData <- rawMeltData[rawMeltData$Name != "Name",]

# Remove samples with no nbr of Tm value
rawMeltData <- rawMeltData[!is.na(rawMeltData$No..Tm.SYBR),]

# Correct fileName to match previous files
rawMeltData$fileName <- gsub("_MeltCurve", ".XLS", rawMeltData$fileName)

# Manual error : on plate 13.06.2018
rawMeltData$Name <- as.character(rawMeltData$Name)
rawMeltData[rawMeltData$fileName %in% "QPCR13.06.2018correct.XLS.csv" &
              rawMeltData$Pos %in% paste0("F", 7:12), "Name"] <- "CEWE_AA_0367"

rawData <- merge(rawData, 
                 rawMeltData[c("fileName", "Name", "Pos", "No..Tm.SYBR", "Tm.x..SYBR")],
                 by = c("fileName", "Name", "Pos"), all = T)

# Remove controls
rawData <- rawData[!rawData$Name %in% c("water", "NTC", "Noiseband","OFF", "Quantification", "SYBR","n/a" ) &
                     !rawData$Pos %in% c("Threshold level", "Baseline setting"),]

# Remove this plate as all samples were redone, better, later + no melt curve
rawData <- rawData[!rawData$fileName %in% "QPCR04.04.2018.XLS.csv",]

# No melt curves for them:
unique(rawData[is.na(rawData$No..Tm.SYBR),"fileName"])

# Remove stupid empty rows
rawData <- rawData[!is.na(rawData$fileName),]

# manual corrections on sample names
rawData$Name[rawData$Name %in% "359"] <- "ILWE_AA_0359"
rawData$Name[rawData$Name %in% "ILWE_AA_242"] <- "ILWE_AA_0242"
rawData$Name[rawData$Name %in% "CEWE_AA_0330#"] <- "CEWE_AA_0330"
rawData$Name[rawData$Name %in% "ILWE_AA_0,48"] <- "ILWE_AA_0348"
rawData$Name[rawData$Name %in% "ILWE_AA_0134"] <- "ILWE_AA_0379"
rawData[rawData$fileName %in% "QPCR07.06.2018part2.XLS.csv" &
          rawData$Pos %in% paste0("C", 4:12),"Name"] <- gsub("0410", "0409",
                                                             rawData[rawData$fileName %in% "QPCR07.06.2018part2.XLS.csv" &
                                                                       rawData$Pos %in% paste0("C", 4:12),"Name"])
rawData[rawData$fileName %in% "QPCR07.06.2018part2.XLS.csv" &
          rawData$Pos %in% paste0("D", 7:12),"Name"] <- "ILWE_AA_0410" 
rawData[rawData$Pos %in% c("B4", "B5", "B6") &
          rawData$fileName %in% "QPCR14.06.2018.XLS.csv", "Name"] <- "CEWE_AA_0424"

# Add full name of sample (tissue + mouseID + eimeriaOrmouse primers + plate)
rawData$fullName <- paste0(rawData$Name, "_", rawData$Target.SYBR, "_",  rawData$fileName)

# Important factor to numbers here:
rawData$Tm.x..SYBR <- as.numeric(as.character(rawData$Tm.x..SYBR))
rawData$Ct.SYBR <- as.numeric(as.character(rawData$Ct.SYBR))
rawData$Ct.Mean.SYBR <- as.numeric(as.character(rawData$Ct.Mean.SYBR))

rawData$triplicate <- paste(rawData$Name, rawData$fileName)

# Add tissue and Mouse_ID
x <- strsplit(as.character(rawData$Name), "_", 1)
rawData$tissue <- sapply( x, "[", 1)
rawData$Mouse_ID <- paste0("AA_", sapply( x, "[", 3))
rm(x)
##### END Prepare data #####

## heatmap to follow
library(ggplot2)

myTiles <- function(df){
  dat_long <- data.frame(tissue_target = paste(rawData$tissue, rawData$Target.SYBR), 
                         variable = rawData$Mouse_ID,
                         value = 0)
  
  dat_long <- unique(dat_long)
  dat_long[paste(dat_long$tissue_target, dat_long$variable) %in% 
             paste(df$tissue, df$Target.SYBR, df$Mouse_ID), "value"] <- 1
  
  # discrete vs continuous
  dat_long$value <- factor(dat_long$value)
  
  gg <- ggplot(dat_long)
  # fill + legend, gray border
  gg <- gg + geom_tile(aes(x = variable, y = tissue_target, fill = value),
                       color="#7f7f7f")
  # custom fill colors
  gg <- gg + scale_fill_manual(values=c("grey", "green"))
  # no labels
  gg <- gg + labs(x=NULL, y=NULL)
  # remove some chart junk
  gg <- gg + theme_bw() + theme(panel.grid=element_blank(),
                                panel.border=element_blank(),
                                axis.text.x = element_text(angle = 45, hjust = 1, size=3) )
  return(list(gg, table(dat_long$tissue_target, dat_long$value)))
}

myTiles(rawData)

##### Clean data #####
## 1. Which samples are complete (mouse+eimeria triplicate, low sd)
# triplicate = same file, same name, same target, same mean

library(dplyr)

sumOKData <- rawData %>% 
  group_by(fullName) %>% 
  summarise(count = n()) %>% 
  data.frame()

OKData <- rawData[!rawData$fullName %in% 
                   sumOKData[sumOKData$count < 3 , ]$fullName, ]
myTiles(OKData)

# 2. Let's have a look at our melting curve information...
ggplot(OKData, aes(Tm.x..SYBR, col = Target.SYBR)) +
  geom_histogram(fill = NA, binwidth = .1) +
  geom_vline(xintercept = 76.3)+
  facet_grid(.~Target.SYBR) +
  theme_bw()

# Plot Eimeria melting curve peak vs Ct to see
ggplot(OKData, aes(x = Tm.x..SYBR, y = Ct.Mean.SYBR)) +
  geom_point() +
  facet_grid(.~Target.SYBR) +
  theme_bw() +
  geom_vline(xintercept = 76, col = "red") +
  ggtitle("Plot after phase 1. keeping all triplicates")

## TO SEE IF WE USE:

# # based on these plots, I would take 10-50 Tm for eimeria positive, and 50 to 80 for mouse positive
# OKData$meltingcurveStatus[OKData$Target.SYBR %in% "eimeria"] <- "NOeimeriaDNA"
# OKData$meltingcurveStatus[
#   as.numeric(as.character(OKData$Tm.x..SYBR)) < 76.3 & 
#     OKData$Target.SYBR %in% "eimeria"] <- "eimeriaDNA"
# OKData$meltingcurveStatus[OKData$Target.SYBR %in% "mouse"] <- "NOmouseDNA"
# OKData$meltingcurveStatus[
#   as.numeric(as.character(OKData$Tm.x..SYBR)) > 76.3 & 
#     OKData$Target.SYBR %in% "mouse"] <- "mouseDNA"
# OKData$meltingcurveStatus[OKData$fileName %in% "QPCR30.05.2018.XLS.csv"] <- "missingInfo"

# 3. how are the sd? keep < 3
badSd <- OKData[!is.na(OKData$Ct.Dev..SYBR) & OKData$Ct.Dev..SYBR > 3 ,]
OKData <- OKData[OKData$Ct.Dev..SYBR <= 3 ,]
myTiles(OKData)

## TO SEE IF WE USE:

# # 4. Remove mice DNA if not present 3 times
# checkOKData <- OKData[OKData$Target.SYBR %in% "mouse",] %>%
#   group_by(fullName) %>%
#   summarize(mouseMeltingCurveStatus = length(table(meltingcurveStatus)),
#             isMousePositive = meltingcurveStatus[1]) %>%
#   data.frame()
# 
# toremove <- checkOKData$fullName[checkOKData$mouseMeltingCurveStatus %in% 2 | 
#   checkOKData$isMousePositive %in% "NOmouseDNA"]
# 
# OKData <- OKData[!OKData$fullName %in% toremove,]

# myTiles(OKData)

## TO SEE IF WE USE:

# 5. Remove also if Emeria DNA is sometimes here, sometimes not. Should be 3 positive or 3 negative!!
# checkOKData <- OKData[OKData$Target.SYBR %in% "eimeria",] %>% 
#   group_by(fullName) %>% 
#   summarize(eimeriaMeltingCurveStatus = length(table(meltingcurveStatus))) %>% 
#   data.frame()
# 
# OKData <- OKData[!OKData$fullName %in% checkOKData$fullName[checkOKData$eimeriaMeltingCurveStatus > 1],]

# myTiles(OKData)

# 4. Calculate delta ct
calculateDeltaCt <- function(df){
  # remove individual values for position and Ct, so that we avoid useless repeats (we care about Ct Sd and Mean)
  
  
  # df <- df[names(df)[!names(df) %in% c("Ct.SYBR", "Pos", "Tm.x..SYBR")]]

  ### To redo later on 
  
  ## Calculate deltaCt per plate
  sumDataMouse <- df[df$Target.SYBR %in% "mouse",]
  sumDataEimeria <- df[df$Target.SYBR %in% "eimeria",]
  mergedData <- merge(sumDataEimeria, sumDataMouse,
                      by = c("Mouse_ID", "fileName", "tissue", "Name"))
  mergedData <- unique(mergedData)
  mergedData$deltaCt <- as.numeric(as.character(mergedData$Ct.Mean.SYBR.y)) - 
    as.numeric(as.character(mergedData$Ct.Mean.SYBR.x)) # DeltaCt EIMERIA minus MOUSE
  return(mergedData)
}

OKData <- calculateDeltaCt(OKData)
BadData <- calculateDeltaCt(badSd)

myTiles2 <- function(df){
  dat_long <- data.frame(tissue = rawData$tissue, 
                         Mouse_ID = rawData$Mouse_ID,
                         value = 0)
  
  dat_long <- unique(dat_long)
  dat_long[paste(dat_long$tissue, dat_long$Mouse_ID) %in% 
             paste(df$tissue, df$Mouse_ID), "value"] <- 1
  
  # discrete vs continuous
  dat_long$value <- factor(dat_long$value)
  
  gg <- ggplot(dat_long)
  # fill + legend, gray border
  gg <- gg + geom_tile(aes(x = Mouse_ID, y = tissue, fill = value),
                       color="#7f7f7f")
  # custom fill colors
  gg <- gg + scale_fill_manual(values=c("grey", "green"))
  # no labels
  gg <- gg + labs(x=NULL, y=NULL)
  # remove some chart junk
  gg <- gg + theme_bw() + theme(panel.grid=element_blank(),
                                panel.border=element_blank(),
                                axis.text.x = element_text(angle = 45, hjust = 1, size=3) )
  return(list(gg, table(dat_long$tissue, dat_long$value)))
}

myTiles2(OKData)
myTiles2(BadData)

OKData$dataset <- "OKData"
BadData$dataset <- "BadData"
AllData <- rbind(OKData, BadData)

# Plot Eimeria melting curve peak vs Ct to see
ggplot(AllData, aes(x = Tm.x..SYBR.y, y = Ct.Mean.SYBR.y)) +
  geom_point(aes(col = dataset), size = 3) +
  theme_bw()

## AND PLOT
# ggplot(OKData, aes(meltingcurveStatus.x, deltaCt)) +
#   geom_boxplot() +
#   geom_jitter() +
#   theme_bw() +
#   geom_hline(yintercept = 6, col = "red")

## And save
OKData$year <- 2017

finalDataClean <- OKData["Mouse_ID"]
# Add CEWE
finalDataClean <- merge(finalDataClean,
                        OKData[OKData$tissue %in% "CEWE", c("Mouse_ID", "deltaCt")],
                           all.x = T)
names(finalDataClean)[names(finalDataClean) %in% "deltaCt"] <- "delta_ct_cewe"

# Add ILWE
finalDataClean <- merge(finalDataClean,
                        OKData[OKData$tissue %in% "ILWE", c("Mouse_ID", "deltaCt")],
                           all.x = T)
names(finalDataClean)[names(finalDataClean) %in% "deltaCt"] <- "delta_ct_ilwe"

# Add observer
finalDataClean$observer_qpcr <- "Lorenzo"

# Write out
write.csv(finalDataClean, "../qPCR_2017.csv", row.names = F)

# ###########
# source("../../../R/functions/addPCRresults.R")
# source("../../../R/functions/addqPCRresults.R")
# source("../../../R/functions/addFlotationResults.R")
#
# myFinal <- addPCRresults(finalData, pathtodata = "../Inventory_contents_all.csv")
# myFinal <- addFlotationResults(myFinal, pathtofinalOO = "../FINALOocysts2015to2017.csv",
#                                pathtolorenzodf = "../Eimeria_oocysts_2015&2017_Lorenzo.csv")$new
#
# myFinal <- addqPCRresults(myFinal, pathtoqPCR2016 = "../qPCR_2016.csv", pathtoqPCR2017 = "../qPCR_2017.csv")
#
# myFinal$year[is.na(myFinal$year)] <- myFinal$year.x[is.na(myFinal$year)]
# myFinal$year[is.na(myFinal$year)] <- myFinal$year.y[is.na(myFinal$year)]
# myFinal$year <- as.factor(myFinal$year)
#
# summary(lm(OPG ~ delta_ct_cewe, myFinal))
#
# ggplot(myFinal, aes(y = myFinal$delta_ct_ilwe, x = OPG+1)) +
#   scale_x_log10() +
#   geom_point(aes(col = year), size = 4) +
#   geom_smooth(method = "lm")
#
# ggplot(myFinal, aes(y = myFinal$delta_ct_cewe, x = OPG+1)) +
#   scale_x_log10() +
#   geom_point(aes(col = year), size = 4) +
#   geom_smooth(method = "lm")
#
# # Plot 1 detection methods compared
# ggplot(myFinal, aes(x = PCRstatus, y = OPG + 1)) +
#   scale_y_log10() +
#   # geom_boxplot()+
#   geom_violin() +
#   geom_jitter(aes(col = year), width = .1, size = 2, alpha = .5) +
#   theme_bw()
#
# # Plot 2 detection methods compared
# ggplot(myFinal, aes(x = qPCRstatus, y = OPG + 1)) +
#   scale_y_log10() +
#   geom_boxplot()+
#   geom_jitter(aes(col = year), width = .1, size = 2, alpha = .5) +
#   theme_bw()
#
# # How to set up a limit of detection for qPCR
# ggplot(myFinal, aes(x = OPG > 0, y = delta_ct_cewe)) +
#   geom_boxplot()+
#   geom_point(aes(col = year), size = 2, alpha = .5) +
#   theme_bw()
#
# ggplot(myFinal, aes(x = myFinal$PCRstatus, y = delta_ct_cewe)) +
#   geom_boxplot()+
#   geom_point(aes(col = year), size = 2, alpha = .5) +
#   theme_bw()
#
# myFinal$FlotOrPcr <- "negative"
# myFinal$FlotOrPcr[myFinal$OPG > 0 | myFinal$PCRstatus %in% "positive"] <- "positive"
#
# ggplot(myFinal, aes(x = myFinal$FlotOrPcr, y = delta_ct_cewe)) +
#   geom_violin()+
#   geom_point(aes(col = year), size = 2, alpha = .5) +
#   geom_hline(yintercept = 6) +
#   theme_bw()
#
# ## Which samples have no qPCR?
# # 2017
# myFinal$Mouse_ID[!is.na(myFinal$delta_ct_cewe)]
# which(duplicated(myFinal$Mouse_ID))
#

