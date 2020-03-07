# fecal ELISas HZ19 CEWE

library(httr)
library(RCurl)
library(Rmisc)
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(drc)
library(data.table)

##### add clean tables ELISA 1
E1_std <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/P3_112019_Eim_CEWE_ELISAs/P3_112019_Eim_CEWE_ELISA1_std.csv"
E1_std <- read.csv(text = getURL(E1_std))

E1_samples <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/P3_112019_Eim_CEWE_ELISAs/P3_112019_Eim_CEWE_ELISA1_samples.csv"
E1_samples <- read.csv(text = getURL(E1_samples))

###### use drc to construct standard curve and pinpointprotein content

model1<-drm(OD~Conc,
            fct=LL.4(names=c("Slope", "Lower", "Upper", "ED50")),
            data=E1_std)
plot(model1)

E1<-ED(model1, E1_samples$OD, type="absolute", display=F)
row.names(E1) <- E1_samples$label

points(y=E1_samples$OD,x=E1[,1],col="lightblue",pch=19,cex=2)
text(y =E1_samples$OD, x = E1[,1], labels=E1_samples$label, data=E1, cex=0.9, font=2)

E1 <- data.frame(E1)
colnames(E1)[1] <- "IFNy"
E1 <- dplyr::select(E1, IFNy)
setDT(E1, keep.rownames = TRUE)[]
colnames(E1)[1] <- "label"
E1 <- merge(E1, E1_samples)
E1$OD <- NULL

write.csv(E1, "./Eimeria_Lab/data/3_recordingTables/P3_112019_Eim_CEWE_ELISAs/EP_112019_Eim_CEWE_ELISA1_complete.csv")
write.csv(E1, "C:/Users/Luke Bednar/Documents/Eimeria_Lab/data/3_recordingTables/P3_112019_Eim_CEWE_ELISAs/P3_112019_Eim_CEWE_ELISA1_complete.csv")

############### load in ELISA2

E2_std <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_P3_Eim_CEWE_ELISA1_std.csv"
E2_std <- read.csv(text = getURL(E2_std))

E2_samples <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/P3_112019_Eim_CEWE_ELISAs/P3_112019_Eim_CEWE_ELISA2_samples.csv"
E2_samples <- read.csv(text = getURL(E2_samples))

###### use drc to construct standard curve and pinpointprotein content

model2<-drm(OD~Conc,
            fct=LL.4(names=c("Slope", "Lower", "Upper", "ED50")),
            data=E2_std)
plot(model2)

E2<-ED(model2, E2_samples$OD, type="absolute", display=F)
row.names(E2) <- E2_samples$label

points(y=E2_samples$OD,x=E2[,1],col="lightblue",pch=19,cex=2)
text(y =E2_samples$OD, x = E2[,1], labels=E2_samples$label, data=E2, cex=0.9, font=2)

E2 <- data.frame(E2)
colnames(E2)[1] <- "IFNy"
E2 <- dplyr::select(E2, IFNy)
setDT(E2, keep.rownames = TRUE)[]
colnames(E2)[1] <- "label"
E2 <- merge(E2, E2_samples)
E2$OD <- NULL

write.csv(E2, "./Eimeria_Lab/data/3_recordingTables/P3_112019_Eim_CEWE_ELISAs/EP_112019_Eim_CEWE_ELISA2_complete.csv")
write.csv(E2, "C:/Users/Luke Bednar/Documents/Eimeria_Lab/data/3_recordingTables/P3_112019_Eim_CEWE_ELISAs/P3_112019_Eim_CEWE_ELISA2_complete.csv")

######################### merge all and write out

E <- rbind(E1, E2)
write.csv(E, "C:/Users/Luke Bednar/Documents/Eimeria_Lab/data/3_recordingTables/P3_112019_Eim_CEWE_ELISAs/P3_112019_Eim_CEWE_ELISAs_complete.csv")
