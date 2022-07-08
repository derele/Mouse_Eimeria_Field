library(tidyverse)

alldata18_19_21 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/CEWE_FECES_infection_intensities")


##how many infected due to fecal samples?
nrow(alldata18_19_21[alldata18_19_21$MC.Eimeria.FEC=="TRUE",])

##how many infected due to tissue sample?

nrow(alldata18_19_21[alldata18_19_21$MC.Eimeria=="TRUE",])

###Plot CqMean Faeces qPCR and delta CT CEWE

class(alldata18_19_21$FEC_Eim_Ct)

alldata18_19_21$FEC_Eim_Ct<-as.numeric(alldata18_19_21$FEC_Eim_Ct)

ggplot(data=alldata18_19_21, mapping = aes(x=FEC_Eim_Ct,y=delta_ct_cewe_MminusE,col=MCs,size=log(OPG+1)))+
  geom_point()+
  scale_size(range=c(0.5,3))+
  geom_smooth(method="lm", size=0.5,se=TRUE)+
  theme(axis.title.x = element_text(size=10))+
  theme(axis.title.y = element_text(size = 10))+
  theme(legend.text=element_text(size=7))+
  theme(legend.title =element_text(size=9))+
  theme(legend.spacing.y = unit(0.05,"cm"))+
  theme(panel.background=element_rect(fill="gray95",color="gray40"))+
  labs(x="Ct from qPCR with faecal DNA",y="deltaCt from qPCR with tissue DNA",color="Infected",size="log(OPG)")+
  guides(color=guide_legend(override.aes = list(size=2)))


ggplot(data=alldata18_19_21,mapping=aes(x=delta_ct_cewe_MminusE,y=FEC_Eim_Ct,col=MCs,size=log(OPG+1)))+
  geom_point()+
  scale_size(range=c(0.5,3))+
  geom_smooth(method="lm", size=0.5,se=TRUE)+
  theme(axis.title.x = element_text(size=10))+
  theme(axis.title.y = element_text(size = 10))+
  theme(legend.text=element_text(size=7))+
  theme(legend.title =element_text(size=9))+
  theme(legend.spacing.y = unit(0.05,"cm"))+
  theme(panel.background=element_rect(fill="gray95",color="gray40"))+
  labs(x="delta Ct from qPCR with tissue DNA",y="Ct from qPCR with faecal DNA",color="Infected",size="log(OPG)")+
  guides(color=guide_legend(override.aes = list(size=2)))

###prediction of amount of oocyst equivalents FEC

Standard_curve<-read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/HZ21_CEWE_EqPCR_SC.csv")

## bringing standard curve into right format
colnames(Standard_curve)[colnames(Standard_curve)%in%"Ct_Mean_Eim"] <-"FEC_Eim_Ct"

colnames(Standard_curve)[colnames(Standard_curve)%in%"Mouse_ID"] <-"Oocysts"

Standard_curve$Oocysts[Standard_curve$Oocysts%in%c("V10_0","V10_1","V10_2","V10_3","V10_4","V10_5","V10_6")]<-c("0","1","2","3","4","5","6")


Standard_curve<- Standard_curve%>%filter(Standard_curve$FEC_Eim_Ct>0)


Eimeria_detection<-alldata18_19_21

##design a linear model
class(Standard_curve$Oocysts)

Standard_curve$Oocysts<-as.numeric(Standard_curve$Oocysts)

Eimeria_detection_FEC<-Eimeria_detection %>% dplyr::select(Mouse_ID,FEC_Eim_Ct)

Linear_model<-lm(Oocysts~FEC_Eim_Ct,data=Standard_curve)

##prediction
class(Eimeria_detection_FEC$FEC_Eim_Ct)
class(Standard_curve$FEC_Eim_Ct)
Eimeria_detection_FEC$FEC_Eim_Ct<-as.numeric(Eimeria_detection_FEC$FEC_Eim_Ct)

Oocysts_Predict_Eimeria_FEC<-10^predict(Linear_model,newdata=Eimeria_detection_FEC)

Eimeria_detection_FEC<-data.frame(Eimeria_detection_FEC,Oocysts_Predict_Eimeria_FEC)

###prediction of amount of oocyst equivalents tissue

Standard_curve2<-read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/HZ21_CEWE_EqPCR_SC.csv")

## bringing standard curve into right format
colnames(Standard_curve2)[colnames(Standard_curve2)%in%"Ct_Mean_Eim"] <-"delta_ct_cewe_MminusE"

colnames(Standard_curve2)[colnames(Standard_curve2)%in%"Mouse_ID"] <-"Oocysts"

Standard_curve2$Oocysts[Standard_curve2$Oocysts%in%c("V10_0","V10_1","V10_2","V10_3","V10_4","V10_5","V10_6")]<-c("0","1","2","3","4","5","6")


Standard_curve2<- Standard_curve2%>%filter(Standard_curve2$delta_ct_cewe_MminusE>0)

##design a linear model

class(Standard_curve2$Oocysts)

Standard_curve2$Oocysts<-as.numeric(Standard_curve2$Oocysts)

Eimeria_detection_CEWE<-Eimeria_detection%>% dplyr::select(Mouse_ID,delta_ct_cewe_MminusE)

Linear_model2<-lm(Oocysts~delta_ct_cewe_MminusE,data=Standard_curve2)


##prediction

Oocysts_Predict_Eimeria_CEWE<-10^predict(Linear_model2,newdata=Eimeria_detection_CEWE)

Eimeria_detection_CEWE<-data.frame(Eimeria_detection_CEWE,Oocysts_Predict_Eimeria_CEWE)

### merging both predictions together

Oocyst_Predict<-merge(Eimeria_detection_FEC,Eimeria_detection_CEWE,by="Mouse_ID",all = TRUE)

sum(duplicated(Oocyst_Predict$Mouse_ID))

Oocyst_Predict<-na.omit(Oocyst_Predict)


##removing not needed dataframes

rm(Eimeria_detection)
rm(Eimeria_detection_CEWE)
rm(Eimeria_detection_FEC)
rm(Standard_curve)
rm(Standard_curve2)

## removing not needed non dataframes

rm(Oocysts_Predict_Eimeria_CEWE)
rm(Oocysts_Predict_Eimeria_FEC)


### plotting predicted oocysts against opgs

## merging Oocyst predict with alldata18_19_21
colnames(Oocyst_Predict)

Oocyst_Predict<-Oocyst_Predict%>% dplyr::select(Mouse_ID,Oocysts_Predict_Eimeria_FEC,Oocysts_Predict_Eimeria_CEWE)

alldata_18_19_21<-merge(alldata18_19_21,Oocyst_Predict,all = TRUE)

colnames(alldata_18_19_21)


alldata_18_19_21<-alldata_18_19_21 %>%
  dplyr::select(Mouse_ID,Year,FEC_Eim_Ct,MC.Eimeria.FEC,Oocysts_Predict_Eimeria_FEC,delta_ct_cewe_MminusE,MC.Eimeria,Oocysts_Predict_Eimeria_CEWE,OPG,MCs)
library(xlsx)

#write.xlsx(alldata_18_19_21,file="alldata_18_19_21.xlsx")

##setting all MC_Eim_FEC FALSE as 0

alldata_18_19_21$Oocysts_Predict_Eimeria_FEC[alldata_18_19_21$MC.Eimeria.FEC=="FALSE"]<-0

##plot

ggplot(data=alldata_18_19_21,mapping=aes(x=log10(OPG),y=log10(Oocysts_Predict_Eimeria_FEC),col=MCs))+
  geom_point()+
  geom_smooth(method="lm", size=0.5,col="black",se=TRUE)+
  theme(axis.title.x = element_text(size=10))+
  theme(axis.title.y = element_text(size = 10))+
  theme(legend.text=element_text(size=7))+
  theme(legend.title =element_text(size=9))+
  theme(legend.spacing.y = unit(0.05,"cm"))+
  theme(legend.background = element_blank())+
  theme(panel.background=element_rect(fill="gray95",color="gray40"))+
  labs(title=NULL,subtitle="a",x="counted oocysts (OPG's)",y="predicted Oocysts from ct")+
  guides(color=guide_legend(override.aes = list(size=2)))

ggplot(data=alldata_18_19_21,mapping=aes(x=log10(Oocysts_Predict_Eimeria_FEC),y=log10(OPG),))+
  geom_point(col="forestgreen")+
  scale_size(range=c(0.5,3))+
  geom_smooth(method="lm", size=0.5,col="black",se=TRUE)+
  theme(axis.title.x = element_text(size=10))+
  theme(axis.title.y = element_text(size = 10))+
  theme(legend.text=element_text(size=7))+
  theme(legend.title =element_text(size=9))+
  theme(legend.spacing.y = unit(0.05,"cm"))+
  theme(legend.background = element_blank())+
  theme(panel.background=element_rect(fill="gray95",color="gray40"))+
  guides(color=guide_legend(override.aes = list(size=2)))

## check which oocyst predict 0 have oocysts

predicted_zero<-alldata_18_19_21%>%select(Mouse_ID,Oocysts_Predict_Eimeria_FEC,OPG,FEC_Eim_Ct)

predicted_zero<-predicted_zero%>%filter(Oocysts_Predict_Eimeria_FEC==0&OPG>0)

predicted_zero<-predicted_zero%>%unique()

Tm1<-FEC_18_19_21%>%select(Mouse_ID,Tm1)

predicted_zero<-merge(Tm1,predicted_zero)

### Check residuals of counted oocysts and predicted oocysts

sum(is.na(alldata_18_19_21$Oocysts_Predict_Eimeria_FEC))

sum(is.na(alldata_18_19_21$OPG))

alldata_18_19_21$OPG[is.na(alldata_18_19_21$OPG)]<-0

alldata_18_19_21$Oocysts_Predict_Eimeria_FEC[is.na(alldata_18_19_21$Oocysts_Predict_Eimeria_FEC)]<-0

rownames(alldata_18_19_21)<-make.unique(alldata_18_19_21$Mouse_ID)

Residual_model<-lm(formula=Oocysts_Predict_Eimeria_FEC~OPG,data=alldata_18_19_21)

Predicted<-predict(Residual_model,newdata=alldata_18_19_21)

RES<-residuals(Residual_model,newdata=alldata_18_19_21)

alldata_18_19_21<-data.frame(alldata_18_19_21,Predicted)

alldata_18_19_21<-merge(alldata_18_19_21,data.frame(RES),by=0,)

plot(RES~Oocysts_Predict_Eimeria_FEC,data=alldata_18_19_21)

ggplot(data=alldata_18_19_21, aes(x=OPG,y=RES))+
  geom_smooth(method="lm",se=F,col="black",size=0.5)+
  geom_point(aes(color=abs(Oocysts_Predict_Eimeria_FEC),size=abs(RES)))+
  scale_size(range=c(1,5))+
  scale_color_continuous(low="indianred2",high="indianred4")+
  labs(title=NULL,subtitle="b",x="Counted oocysts(OPG's)",y="Residuals",color="Predicted genome equivaltents faeces",size="Residuals")+
  guides(color=guide_legend(override.aes=list(size=3)))+
  theme(legend.text=element_text(size=7))+
  theme(legend.title =element_text(size=9))+
  theme(panel.background=element_rect(fill="gray95",color="gray40"))+
  theme(legend.spacing.y = unit(0.05,"cm"))

ggplot(data=alldata_18_19_21, aes(x=OPG,y=Oocysts_Predict_Eimeria_FEC))+
  geom_smooth(method="lm",se=F,col="black",size=0.5)+
  geom_point(aes(color=abs(Oocysts_Predict_Eimeria_CEWE),size=abs(RES)))+
  scale_size(range=c(1,5))+
  scale_color_continuous(low="indianred2",high="indianred4")+
  labs(x="Counted oocysts (OPG's)",y="Predicted genome equivalents faeces",color="Predicted genome equivaltents tissue",size="Residuals")+
  guides(color=guide_legend(override.aes=list(size=3)))+
  theme(legend.text=element_text(size=7))+
  theme(legend.title =element_text(size=9))+
  theme(panel.background=element_rect(fill="gray95",color="gray40"))+
  theme(legend.spacing.y = unit(0.05,"cm"))  


ggplot(data=alldata_18_19_21, aes(x=Oocysts_Predict_Eimeria_FEC,y=RES))+
  geom_smooth(method="lm",se=F,col="black",size=0.5)+
  geom_point(aes(color=abs(Oocysts_Predict_Eimeria_CEWE),size=abs(RES)))+
  scale_size(range=c(1,5))+
  scale_color_continuous(low="indianred2",high="indianred4")+
  guides(color=guide_legend(override.aes=list(size=3)))+
  theme(legend.text=element_text(size=7))+
  theme(legend.title =element_text(size=9))+
  theme(panel.background=element_rect(fill="gray95",color="gray40"))+
  labs(x="Predicted genome equivalents faeces",y="Residuals",color="Predicted genome equivaltents tissue",size="Residuals")+
  theme(legend.spacing.y = unit(0.05,"cm"))

