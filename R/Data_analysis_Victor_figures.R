##Data Analysis and figures 
##Venn-diagram colon content sequences
##Making a data base with all samples 2014 - 2017

# load packages
library("VennDiagram")
library(grid)
library(gridExtra)
library(ggmap)
library(ggrepel)
library(RCurl)
library(ggmap)
library(pegas)
library(msa)
library(RColorBrewer)
library(Biostrings)
library(IRanges)
library(XVector)

# load HMHZ functions
source("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/R/HMHZ_Functions.R")

# Local instruction, only for Victor's computer ;)
setwd("/home/victor/Dokumente/Git_projects/Mouse_Eimeria_Databasing/raw_data/Eimeria_detection/")

pcr.data <- read.csv("Inventory_contents_all.csv")

pcr.data <- pcr.data[-c(5:11,14, 16, 18, 20,21)]

colnames(pcr.data) <- c("Year", "Transect", "Code", "Mouse_ID", "Flot", "Ap5", "n18S_Seq", 
                        "COI_Seq", "ORF470_Seq")

# Data exploration
length(which(pcr.data$Ap5 == "positive" | pcr.data$Flot == "positive"))
length(which(pcr.data$Ap5 == "positive"))
length(which(pcr.data$Flot == "positive"))
length(which(pcr.data$n18S_Seq == "positive"))
length(which(pcr.data$COI_Seq == "positive"))
length(which(pcr.data$ORF470_Seq == "positive"))

# Venn diagram
grid.newpage()
draw.triple.venn(area1 = nrow(subset(pcr.data, n18S_Seq == "positive")), area2 = nrow(subset(pcr.data, COI_Seq == "positive")), 
                 area3 = nrow(subset(pcr.data, ORF470_Seq == "positive")), n12 = nrow(subset(pcr.data, n18S_Seq == "positive"&COI_Seq== "positive")), 
                 n23 = nrow(subset(pcr.data, COI_Seq== "positive"&ORF470_Seq =="positive")), n13 = nrow(subset(pcr.data, n18S_Seq == "positive"&ORF470_Seq == "positive")), 
                 n123 = nrow(subset(pcr.data, n18S_Seq == "positive"&ORF470_Seq == "positive"&COI_Seq=="positive")), category = c("18S", "COI", "ORF470"), 
                 lty = rep(1,3), col = c("dodgerblue4", "firebrick3", "darkgreen"), lwd = rep(2,3),
                 fill = c("dodgerblue4", "firebrick3", "darkgreen"), alpha = c(0.3, 0.3, 0.3), cex =2, cat.cex = 2.5, cat.default.pos = 'outer', 
                 cat.col = c("dodgerblue4", "firebrick3", "darkgreen"))

##Venn-diagram tissue
## ONLY local Victor so far
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
filt.data <- rawdata[names(rawdata) %in% c("Mouse_ID", "Transect", "Sex", "Longitude", "Latitude", "Year", "HI")]

##Set working directory LOCAL
setwd("/home/victor/Dokumente/Sequences/Samples_2017")

# Load haplogroups
haplo <- read.csv("Haplogroups_14_17.csv")
haplo <- haplo[-c(1:3,5:10)]
colnames(haplo) <- c("Mouse_ID", "Haplogroup")
haplo$Mouse_ID <- gsub(pattern = " ", replacement = "", x = haplo$Mouse_ID)

# Replace extra spaces
filt.data$Mouse_ID <- gsub(pattern = " ", replacement = "", x = filt.data$Mouse_ID)
filt.data$Mouse_ID <- gsub(pattern = "Sk3173", replacement = "SK_3173", x = filt.data$Mouse_ID)

# Merge Data bases
Total.DB <- merge(filt.data, haplo , by="Mouse_ID")

Haplo.map(Total.DB, margin = .3)

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

# Keep good data
DB_Results <- DB_Results[!is.na(DB_Results$Longitude) & !is.na(DB_Results$Latitude) & DB_Results$Longitude < 16,]
DB_Results <- DB_Results[complete.cases(DB_Results$Haplogroup),]

##Stablish the area 
Haplo.map(df = DB_Results,  margin = .2)

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

## Written on the 30th August 2017 by Alice Balard/Victor Hugo Jarquin

#####################################################################################

source("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/R/HMHZ_Functions.R")

#let's cheat, I have to git that back after:

# Function creating haplotypes:

##With ORF
VicAlign <- "/home/victor/Dokumente/Sequences/Samples_2014_2015/Haplotype_network/ORF470_Haplotype_300817.fasta"
VicAlign2 <- "/home/victor/Dokumente/Sequences/Samples_2014_2015/Haplotype_network/ORF470_Haplotype_300817.2.fasta"
## ORF (E. falciformis/E. vermiformis/E. separata as references)
VicAlign6 <- "/home/victor/Dokumente/Sequences/Samples_2014_2015/Haplotype_network/ORF470_Haplotype_310817.fasta"

##With 18S
VicAlign3 <- "/home/victor/Dokumente/Sequences/Samples_2014_2015/Haplotype_network/18S_Haplotype_300817.2.fasta"
## 18S (E. falciformis/E. vermiformis/E. ferrisi/E. separata as references)
VicAlign7 <- "/home/victor/Dokumente/Sequences/Samples_2014_2015/Haplotype_network/n18S_Haplotype_310817.2.fasta"

##With COI
VicAlign4 <- "/home/victor/Dokumente/Sequences/Samples_2014_2015/Haplotype_network/COI_Haplotype_300817.2.fasta"
VicAlign8 <- "/home/victor/Dokumente/Sequences/Samples_2014_2015/Haplotype_network/COI_Haplotype_310817.fasta"

##Concatenated Haplotype
VicAlign5 <- "/home/victor/Dokumente/Sequences/Samples_2014_2015/Haplotype_network/Concatenated_haplotype_300817.fasta"
VicAlign9 <- "/home/victor/Dokumente/Sequences/Samples_2014_2015/Haplotype_network/Concatenated_Haplotype_310817.fasta"
VicAlign10 <- "/home/victor/Dokumente/Sequences/Samples_2014_2015/Haplotype_network/Concatenated_Haplotype_060917.fasta"
VicAlign11 <- "/home/victor/Dokumente/Sequences/Samples_2014_2015/Haplotype_network/Concatenated_Haplotype_060917.2.fasta"

##ORF
d <- ape::read.dna(VicAlign, format = "fasta")
d <- ape::read.dna(VicAlign2, format = "fasta")
d <- ape::read.dna(VicAlign6, format = "fasta")
##18S  
d <-ape::read.dna(VicAlign3, format = "fasta")
d <-ape::read.dna(VicAlign7, format = "fasta")
##COI
d <- ape::read.dna(VicAlign4, format = "fasta")
d <- ape::read.dna(VicAlign8, format = "fasta")

##Concatenated
d <- ape::read.dna(VicAlign5, format = "fasta")
d <- ape::read.dna(VicAlign9, format = "fasta")
d <- ape::read.dna(VicAlign10, format = "fasta")
d <- ape::read.dna(VicAlign11, format = "fasta")

## Turn the row names into corresponding HI :
seqnames <- data.frame(Mouse_ID = labels(d))
seqnames
matchname <- merge(seqnames, filt.data, by = "Mouse_ID", all.x = TRUE, sort = FALSE)
matchname
matchname$HI <- round(matchname$HI, 2)
matchname
# Give ref names for ref sequences
matchname[grep("E_", matchname$Mouse_ID), ]$HI <- as.character(matchname[grep("E_", matchname$Mouse_ID), ]$Mouse_ID)
# Give name for NA data
matchname[which(is.na(matchname$HI)),]$HI <- "unknow_yet"
matchname
# Merge all
match <- merge(seqnames, matchname, by = "Mouse_ID", sort = FALSE)
match
## Change name of the sequence by corresponding HI :
rownames(d) <- match$HI

e <- dist.dna(d)

h <- pegas::haplotype(d)

h <- sort(h, what = "label")

(net <- pegas::haploNet(h))

ind.hap<-with(
  stack(setNames(attr(h, "index"), rownames(h))),
  table(hap=ind, pop=rownames(d)[values])
)

## Set colors:
##Change the number of colors acording to the number of ref sequences + Unknown 
mycols <- c(colorRampPalette(c("blue", "red"))(ncol(ind.hap) - 6),  "green", "darkgreen", "darkorange", "yellow", "hotpink",  "grey")

plot(net, size = attr(net, "freq"), pie = ind.hap, fast = TRUE, scale.ratio = 1, 
     cex = 1, bg = mycols) 
legend(15, -10, colnames(ind.hap), fill = mycols, pch=3, ncol=2, cex = 0.6) #+ scale_color_gradient("Hybrid\nindex", high="red",low="blue")
#}
