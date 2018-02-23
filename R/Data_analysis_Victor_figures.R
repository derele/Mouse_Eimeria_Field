##Data Analysis and figures 

##Venn-diagram colon content sequences

##Making a data base with all samples 2014 - 2016

setwd("/home/victor/Dokumente/Git_projects/Mouse_Eimeria_Databasing/raw_data/Eimeria_detection/")


pcr.data <- read.csv("Inventory_contents_all.csv")

pcr.data <- pcr.data[-c(5:11,14, 16, 18, 20,21)]


colnames(pcr.data) <- c("Year", "Transect", "Code", "Mouse_ID", "Flot", "Ap5", "n18S_Seq", 
                        "COI_Seq", "ORF470_Seq")


length(which(pcr.data$Ap5 == "positive" | pcr.data$Flot == "positive"))
length(which(pcr.data$Ap5 == "positive"))
length(which(pcr.data$Flot == "positive"))
length(which(pcr.data$n18S_Seq == "positive"))
length(which(pcr.data$COI_Seq == "positive"))
length(which(pcr.data$ORF470_Seq == "positive"))

library("VennDiagram")
library(grid)
library(gridExtra)

grid.newpage()
draw.triple.venn(area1 = nrow(subset(pcr.data, n18S_Seq == "positive")), area2 = nrow(subset(pcr.data, COI_Seq == "positive")), 
                 area3 = nrow(subset(pcr.data, ORF470_Seq == "positive")), n12 = nrow(subset(pcr.data, n18S_Seq == "positive"&COI_Seq== "positive")), 
                 n23 = nrow(subset(pcr.data, COI_Seq== "positive"&ORF470_Seq =="positive")), n13 = nrow(subset(pcr.data, n18S_Seq == "positive"&ORF470_Seq == "positive")), 
                 n123 = nrow(subset(pcr.data, n18S_Seq == "positive"&ORF470_Seq == "positive"&COI_Seq=="positive")), category = c("18S", "COI", "ORF470"), 
                 lty = rep(1,3), col = c("dodgerblue4", "firebrick3", "darkgreen"), lwd = rep(2,3),
                 fill = c("dodgerblue4", "firebrick3", "darkgreen"), alpha = c(0.3, 0.3, 0.3), cex =2, cat.cex = 2.5, cat.default.pos = 'outer', 
                 cat.col = c("dodgerblue4", "firebrick3", "darkgreen"))



##Venn-diagram tissue

setwd("/home/victor/Dokumente/Sequences/Tissue/")

tissue.data <- read.csv("Haplotypes_Tissue_compleat.csv")

colnames(tissue.data) <- c("Mouse_ID", "Tissue", "n18S_Seq", "n18S_Haplo",
                        "COI_Seq", "COI_Haplo", "ORF470_Seq", "ORF470_Haplo", "Haplogroup")

length(which(tissue.data$n18S_Seq == "positive"))
length(which(tissue.data$COI_Seq == "positive"))
length(which(tissue.data$ORF470_Seq == "positive"))

grid.newpage()
draw.triple.venn(area1 = nrow(subset(tissue.data, n18S_Seq == "positive")), area2 = nrow(subset(tissue.data, COI_Seq == "positive")), 
                 area3 = nrow(subset(tissue.data, ORF470_Seq == "positive")), n12 = nrow(subset(tissue.data, n18S_Seq == "positive"&COI_Seq== "positive")), 
                 n23 = nrow(subset(tissue.data, COI_Seq== "positive"&ORF470_Seq =="positive")), n13 = nrow(subset(tissue.data, n18S_Seq == "positive"&ORF470_Seq == "positive")), 
                 n123 = nrow(subset(tissue.data, n18S_Seq == "positive"&ORF470_Seq == "positive"&COI_Seq=="positive")), category = c("18S", "COI", "ORF470"), 
                 lty = rep(1,3), col = c("dodgerblue4", "firebrick3", "darkgreen"), lwd = rep(2,3),
                 fill = c("dodgerblue4", "firebrick3", "darkgreen"), alpha = c(0.3, 0.3, 0.3), cex =2, cat.cex = 2.5, cat.default.pos = 'outer', 
                 cat.col = c("dodgerblue4", "firebrick3", "darkgreen"))


##Map for all samples (Alice super clean table mice without HI are not included)
rawdata <- read.csv ("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/raw_data/MiceTable_2014to2017.csv")
filt.data <- rawdata[-c(8:28,31:100)]
filt.data <- filt.data[-c(10:20)]

##Set working directory 
setwd("/home/victor/Dokumente/Sequences/Samples_2017")

haplo <- read.csv("Haplogroups_14_17.csv")

haplo <- haplo[-c(1:3,5:10)]

colnames(haplo) <- c("Mouse_ID", "Haplogroup")

##find extra spaces in the lables
tail(as.character(filt.data$Mouse_ID))
tail(as.character(haplo$Mouse_ID))

##Replace extra spaces
filt.data$Mouse_ID <- gsub(pattern = " ", replacement = "", x = filt.data$Mouse_ID)
filt.data$Mouse_ID <- gsub(pattern = "Sk3173", replacement = "SK_3173", x = filt.data$Mouse_ID)

haplo$Mouse_ID <- gsub(pattern = " ", replacement = "", x = haplo$Mouse_ID)

#Merge Data bases

Total.DB <- merge(filt.data, haplo , by="Mouse_ID")

##Stablish the area 
library(ggmap)
area <- get_map(location =
                  c(min(Total.DB$Longitude-0.3),
                    min(Total.DB$Latitude-0.3),
                    max(Total.DB$Longitude+0.3),
                    max(Total.DB$Latitude+0.3)),
                source = "stamen", maptype="toner-lite")


ggmap(area) +
  geom_point(data = Total.DB,
             aes(x = Longitude, y = Latitude, colour = Haplogroup), size = 3,  alpha = 0.7 ) + 
  scale_color_manual(breaks = c("A", "B", "C", "CA", "CB"),
                     values=c("#FF3300","#66CC00", "#FFCC00", "grey", "black"))


##Map for all samples (Victor manually created table all mice included)

##Set working directory 
setwd("/home/victor/Dokumente/Sequences/")

## Get databases
DB_Results <- read.csv("GLM_data.csv")

##find extra spaces in the lables
tail(as.character(DB_Results$Mouse_ID))

##Replace extra spaces
DB_Results$Mouse_ID <- gsub(pattern = " ", replacement = "", x = DB_Results$Mouse_ID)

##Merge with haplogroup by mouse ID

DB_Results <- merge(DB_Results, haplo , by="Mouse_ID")

##Stablish the area 
library(ggmap)
area <- get_map(location =
                  c(min(DB_Results$Longitude-0.3),
                    min(DB_Results$Latitude-0.3),
                    max(DB_Results$Longitude+0.3),
                    max(DB_Results$Latitude+0.3)),
                source = "stamen", maptype="toner-lite")


ggmap(area) +
  geom_point(data = DB_Results,
             aes(x = Longitude, y = Latitude, colour = Haplogroup), size = 3,  alpha = 0.7 ) + 
  scale_color_manual(breaks = c("A", "B", "C", "CA", "CB"),
                     values=c("#FF3300","#66CC00", "#FFCC00", "grey", "black"))


##Violin plot

library("ggplot2")

ggplot(na.omit(Total.DB), aes(x = Haplogroup, y = HI)) + 
  geom_violin(aes(fill = Haplogroup), alpha = 0.7) + 
  scale_fill_manual(values = c("#FF3300","#66CC00", "#FFCC00"))+ geom_jitter(width = 0.1) + 
  labs(title="Distribution of Eimeria haplogroups according to host genotype",x="Haplogroup", y = "HI") + theme_bw() + 
  theme(plot.title = element_text(size=24, face = "bold"), axis.text= element_text(size = 16), axis.title=element_text(size=18,face="bold"))

ggplot(na.omit(merge.data), aes(x = Transect, y = HI)) + 
  geom_violin(aes(fill = Transect), alpha = 0.7) + 
  scale_fill_manual(values = c("#FF3300","#66CC00", "#FFCC00"))+ geom_jitter(width = 0.1) + 
  labs(title="Distribution of HI according to transect",x="Transect", y = "HI") + theme_bw()

##Haplotype network 

