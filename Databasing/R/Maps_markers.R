# MAPS for presentation 17 June 2016 Alice
#setwd("/home/alice/Dokumente/Presentations_Alice/3- 17June2016/")
tabMarkersStatus <- read.csv("/home/alice/Dokumente/Presentations_Alice/3- 17June2016/data2014mouse_vs_eimeria.csv")
tabGenotypes2014 <- read.csv("/home/alice/Dropbox/Eimeria_Wild/mouse_data/HZ14_Mice 31-12-14_genotypes.csv")

# Adapted quickly from Emanuel, to git ASAP!!
library(ggmap)
library(png)
library(grid)

tabGenotypes2014$sample <- paste(tabGenotypes2014$ID, tabGenotypes2014$PIN) # unifromise the sample names in the 2 dataframes
tabGenotypes2014 <- tabGenotypes2014[-1,] # 1st line useless

# Keep only the lines of tabGenotypes2015 found also in tabMarkersStatus
choiceoflines <- tabGenotypes2014$sample %in% tabMarkersStatus$sample
index <- which(choiceoflines) #give the lines where tabGenotypes2015$sample has a status

tabGenotypes2014withStatus <- tabGenotypes2014[index,]

DataGenoStatus <- merge(tabGenotypes2014withStatus,tabMarkersStatus,by = "sample")
DataGenoStatus <- DataGenoStatus[,c("Code","sample","X332","X347","X65","Tsx","Btk","Syap1","X18s","COX_seq","Ap5_seq","type")]

####### Add lat and long to the Dataframe
tabLocal2014 <- read.csv("/home/alice/Dropbox/Eimeria_Wild/mouse_data/HZ14_Mice 31-12-14_localities.csv")
tabLocal2014 <- tabLocal2014 [,c("Latitude","Longitude","Code")]


Datatot <- merge(tabLocal2014, DataGenoStatus, by = "Code", all = TRUE, sort = FALSE)
Datatot <- na.omit(Datatot)

#add size (number of repeat for each place) to Datatot
Number<-as.data.frame(table(Datatot$Code))
Number <-subset(Number,Number$Freq!=0)
names(Number) <- c("Code", "Sample_size")
Datatot <- merge(Datatot, Number, by = "Code", all = TRUE, sort = FALSE)
###########

### Express lat and long better (code Emanuel)

convert.deg <-function(c){
  z <- lapply(strsplit(as.character(c),
                       "° *|' *|(\" *|$)"), as.numeric)
  ## fill to seconds if absent
  zz <- lapply(z, function(X) {
    c(X, rep(0, times = 3 - length(X)))
  })
  
  dec <- lapply(zz, function (x) x[1] + x[2]/60 + x[3]/3600)
  return(unlist(dec))
} 

Datatot$Latitude <- convert.deg(Datatot$Latitude)
Datatot$Longitude <- convert.deg(Datatot$Longitude)

## Could be nice to calculate HI here later

area <- get_map(maptype = "toner",source = "stamen", location = c(min(Datatot$Longitude)-0.2, min(Datatot$Latitude)-0.2, 
                             max(Datatot$Longitude)+0.2,max(Datatot$Latitude)+0.2)) # Get the same map for all plots

# Plot the sampling size:
ggmap(area)+
  geom_point(aes(x = Longitude, y = Latitude, size = Sample_size),
             data = Datatot, alpha = 0.8, color="red")+
  theme(legend.position=c(.1, .9))

# Replace 0 by NA for the markers COI, AP, 18
Datatot$COX_seq[Datatot$COX_seq == 0] <- NA
Datatot$Ap5_seq[Datatot$Ap5_seq == 0] <- NA
Datatot$X18s[Datatot$X18s == 0] <- NA

# Plot
plot18S <-ggmap(area,extent = "device")+
  geom_point(aes(x = Longitude, y = Latitude, color= X18s),
             data = Datatot, alpha = 0.7, size=5, position = position_jitter(0.06,0.06))+
  scale_color_manual(values=c("blue","orange"))+
  theme(legend.position="none")+
  ggtitle("18s")
  
plotCox <- ggmap(area,extent = "device")+
  geom_point(aes(x = Longitude, y = Latitude, color= COX_seq),
             data = Datatot, alpha = 0.7, size=5, position = position_jitter(0.06,0.06))+
  scale_color_manual(values=c("gold","purple"))+
  theme(legend.position="none")+
  ggtitle("COI")

plotAp <-ggmap(area,extent = "device")+
  geom_point(aes(x = Longitude, y = Latitude, color= Ap5_seq),
             data = Datatot, alpha = 0.7, size=5, position = position_jitter(0.06,0.06))+
  scale_color_manual(values=c("green","red"))+
  theme(legend.position="none")+
  ggtitle("Ap5")

require(gridExtra)
grid.arrange(plot18S, plotCox, plotAp, ncol=3)
g <- arrangeGrob(plot18S, plotCox, plotAp, ncol=3)

ggsave("3markers.pdf",g, width = 49, height = 49)
