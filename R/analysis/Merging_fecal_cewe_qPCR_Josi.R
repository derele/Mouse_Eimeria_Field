library(tidyverse)


#### Taken from Josi's Bachelor thesis: 
###FAECAL SAMPLES

#Reading the data product of cleaned data from the fecal qpcr
FEC_18_19_21 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/FEC_EqPCR_DataProduct.csv")

##selecting right TM 

ggplot(data=FEC_18_19_21,aes(x=Tm1,y=stat(density)))+
  geom_histogram(binwidth=0.1,fill="black")+
  geom_density(col="red")+
  geom_vline(xintercept =73.5,col="darkolivegreen3",size=0.5)+
  geom_vline(xintercept=76,col="darkolivegreen3",size=0.5)+
  labs(title= NULL,subtitle="a",x="Tm1 in °C",y="Density")+
  theme(plot.title=element_text(hjust=0.0))+
  theme(axis.text.x=element_text(size=8))+
  theme(panel.background=element_rect(fill="gray95",color="gray40"))+
  scale_x_continuous(breaks=seq(62,93,1.0))

##specifying infected samples

FEC_18_19_21$Tmcorrect <- if_else(between(FEC_18_19_21$Tm1,73.5,76), TRUE, FALSE)

FEC_18_19_21 <- FEC_18_19_21 %>%
  group_by(Mouse_ID)%>%
  mutate(MC.Eimeria.FEC = sum(Tmcorrect)>=3)

##removing doubled rows

qPCR_mouse <- FEC_18_19_21 %>%
  dplyr::select(c(Mouse_ID, FEC_Eim_Ct, Cq.SD, MC.Eimeria.FEC)) %>%
  unique()

qPCR_mouse <- dplyr::rename(qPCR_mouse, Cq.SD.FEC = Cq.SD)

##merging with SOTA

SOTA <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/SOTA_Data_Product.csv")


allData <- merge(SOTA, qPCR_mouse, by="Mouse_ID")


colnames(allData)

###CEWE data

CEWE_18_19_21<-read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/CEWE_EqPCR_Data_Product.csv")


##deciding right tm

ggplot(data=CEWE_18_19_21,aes(x=Tm1_Eim,y=stat(density)))+
  geom_histogram(binwidth=0.05,fill="black")+
  geom_density(col="red")+
  geom_vline(xintercept =73.5,col="darkolivegreen3",size=0.5)+
  geom_vline(xintercept=76,col="darkolivegreen3",size=0.5)+
  labs(title=NULL, subtitle = "b",x="Tm1 in °C", y="Density")+
  theme(plot.title=element_text(hjust=0.0))+
  theme(axis.text.x=element_text(size=8))+
  theme(panel.background=element_rect(fill="gray95",color="gray40"))+
  scale_x_continuous(breaks=seq(72,100,0.5))

##specifying infected samples

CEWE_18_19_21$Tmcorrect <- if_else(between(CEWE_18_19_21$Tm1_Eim,73.5,76), TRUE, FALSE)

CEWE_18_19_21 <- CEWE_18_19_21 %>%
  group_by(Mouse_ID)%>%
  mutate(MC.Eimeria = sum(Tmcorrect)>=3)

##removing doubled rows
CEWE_18_19_21 <- CEWE_18_19_21 %>%
  dplyr::rename(Ct.Eimeria = Ct_Mean_Eim)

qPCR_mouse_CEWE <- CEWE_18_19_21 %>%
  dplyr::select(c(Mouse_ID, Ct.Eimeria, MC.Eimeria)) %>%
  unique()

qPCR_mouse_CEWE[!duplicated(qPCR_mouse_CEWE$Mouse_ID),]->qPCR_mouse_CEWE

colnames(qPCR_mouse_CEWE)
##merging with SOTA

alldata_CEWE <- merge(allData, qPCR_mouse_CEWE, by="Mouse_ID", all = TRUE)

alldata_CEWE<- alldata_CEWE %>% 
  mutate(MC.Eimeria=coalesce(MC.Eimeria.x,MC.Eimeria.y))

alldata_CEWE<- alldata_CEWE %>%
  mutate(Ct.Eimeria=coalesce(Ct.Eimeria.x,Ct.Eimeria.y))

colnames(alldata_CEWE)
alldata_CEWE<-alldata_CEWE %>% 
  dplyr::select(Mouse_ID,Year,MC.Eimeria,delta_ct_cewe_MminusE,FEC_Eim_Ct,OPG,MC.Eimeria.FEC)


### sorting out interesting years 2018,2019,2021

colnames(alldata_CEWE)



alldata18_19_21<- alldata_CEWE %>%
  dplyr::select(Mouse_ID,Year,FEC_Eim_Ct,MC.Eimeria.FEC,delta_ct_cewe_MminusE,MC.Eimeria,OPG)



##paste both MC into one column
colnames(alldata18_19_21)

alldata18_19_21<-alldata18_19_21 %>% 
  dplyr::mutate(MCs=case_when(MC.Eimeria.FEC==TRUE& MC.Eimeria==TRUE~'Fec:TRUE Tissue:TRUE',
                                                        MC.Eimeria.FEC==TRUE& MC.Eimeria==FALSE~'Fec:TRUE Tissue:FALSE',
                                                        MC.Eimeria.FEC==TRUE& is.na(MC.Eimeria)~'Fec:TRUE Tissue:NA',
                                                        MC.Eimeria.FEC==FALSE& MC.Eimeria==TRUE~'Fec:FALSE Tissue:TRUE',
                                                        MC.Eimeria.FEC==FALSE& MC.Eimeria==FALSE~'Fec:FALSE Tissue:FALSE',
                                                        MC.Eimeria.FEC==FALSE& is.na(MC.Eimeria)~'Fec:FALSE Tissue:NA',
                                                        is.na(MC.Eimeria.FEC)& MC.Eimeria==TRUE~'Fec:NA Tissue:TRUE',
                                                        is.na(MC.Eimeria.FEC)&MC.Eimeria==FALSE~'Fec:NA Tissue:FALSE',
                                                        is.na(MC.Eimeria.FEC)&is.na(MC.Eimeria)~"Fec:NA Tissue:NA"))

##removing doubled IDs by combining rows

alldata18_19_21[!duplicated(alldata18_19_21$Mouse_ID),] -> alldata18_19_21

alldata18_19_21$FEC_Eim_Ct[alldata18_19_21$Mouse_ID%in%"AA_0538"]<-"35.35597"

alldata18_19_21$FEC_Eim_Ct[alldata18_19_21$Mouse_ID%in%"AA_0814"]<-"27.57534"



write.csv(alldata18_19_21, "data_input/CEWE_FECES_infection_intensities", row.names=FALSE)


