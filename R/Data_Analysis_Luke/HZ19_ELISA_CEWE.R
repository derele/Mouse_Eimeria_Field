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
HZ_std1 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/ELISAs/HZ19_CEWE_ELISA1_std.csv"))

HZ19_samples1 <- read.csv(text = getURL("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/ELISAs/HZ19_CEWE_ELISA1_samples.csv"))

###### use drc to construct standard curve and pinpointprotein content

model1<-drm(OD~Conc,
            fct=LL.4(names=c("Slope", "Lower", "Upper", "ED50")),
            data=HZ_std1)
plot(model1)

HZ1<-ED(model1, HZ19_samples1$OD, type="absolute", display=F)
row.names(HZ1) <- HZ19_samples1$Mouse_ID

points(y=HZ19_samples1$OD,x=HZ19_samples1[,1],col="lightblue",pch=19,cex=2)
text(y =HZ19_samples1$OD, x = HZ1[,1], labels=HZ19_samples1$Mouse_ID, data=HZ1, cex=0.9, font=2)

HZ1 <- data.frame(HZ1)
colnames(HZ1)[1] <- "IFNy"
HZ1 <- dplyr::select(HZ1, IFNy)
setDT(HZ1, keep.rownames = TRUE)[]
colnames(HZ1)[1] <- "Mouse_ID"
HZ1<- merge(HZ1, HZ19_samples1)
HZ1$OD <- NULL

# write.csv(HZ1, "./Documents/Mouse_Eimeria_Databasing/data/ELISAs/HZ19_CEWE_ELISA1_complete.csv")
write.csv(E1, "C:/Users/Luke Bednar/Mouse_Eimeria_Databasing/data/ELISAs/HZ19_CEWE_ELISA1_complete.csv")

############### load in ELISA2

E2_std <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/ELISAs/HZ19_CEWE_ELISA2_std.csv"
E2_std <- read.csv(text = getURL(E2_std))

E2_samples <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/ELISAs/HZ19_CEWE_ELISA2_samples.csv"
E2_samples <- read.csv(text = getURL(E2_samples))

###### use drc to construct standard curve and pinpointprotein content

model2<-drm(OD~Conc,
            fct=LL.4(names=c("Slope", "Lower", "Upper", "ED50")),
            data=E2_std)
plot(model2)

E2<-ED(model2, E2_samples$OD, type="absolute", display=F)
row.names(E2) <- E2_samples$Mouse_ID

points(y=E2_samples$OD,x=E2[,1],col="lightblue",pch=19,cex=2)
text(y =E2_samples$OD, x = E2[,1], labels=E2_samples$Mouse_ID, data=E2, cex=0.9, font=2)

E2 <- data.frame(E2)
colnames(E2)[1] <- "IFNy"
E2 <- dplyr::select(E2, IFNy)
setDT(E2, keep.rownames = TRUE)[]
colnames(E2)[1] <- "Mouse_ID"


write.csv(E2, "C:/Users/Luke Bednar/Mouse_Eimeria_Databasing/data/ELISAs/HZ19_CEWE_ELISA2_complete.csv")
# write.csv(E2, "C:/Users/Luke Bednar/Documents/Eimeria_Lab/data/3_recordingTables/P3_112019_Eim_CEWE_ELISAs/P3_112019_Eim_CEWE_ELISA2_complete.csv")

######################### merge all and write out

E <- rbind(HZ1, E2)
write.csv(E, "C:/Users/Luke Bednar/Mouse_Eimeria_Databasing/data/ELISAs/HZ19_CEWE_ELISAs_complete.csv")

#################
