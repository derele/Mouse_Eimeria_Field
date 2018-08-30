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
library("ggplot2")
library("BarcodingR")
library(dplyr)
library("ggpubr")
library(multcomp)
library("Publish") 
library("AICcmodavg")
library(Rmisc)
library("ggpubr")

# load HMHZ functions
source("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/R/functions/HMHZ_Functions.R")
source("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/R/functions/makeMiceTable.R")
# Local instruction, only for Victor's computer ;)
setwd("/home/victor/Dokumente/Git_projects/Mouse_Eimeria_Databasing/raw_data/Eimeria_detection/")

#######################Data exploration##############
pcr.data <- read.csv("Inventory_contents_all.csv")

pcr.data <- pcr.data[-c(5:11,14, 16, 18, 20, 21)]

colnames(pcr.data) <- c("Year", "Transect", "Code", "Mouse_ID", "Flot", "Ap5", "n18S_Seq", 
                        "COI_Seq", "ORF470_Seq","Cecum_Dx", "Ileum_Dx")


length(which(pcr.data$Ap5 == "positive" | pcr.data$Flot == "positive"))
length(which(pcr.data$Ap5 %in% "negative"))

length(which(pcr.data$Ap5 == "positive" | pcr.data$Flot == "positive" | pcr.data$Cecum_Dx == "positive" | pcr.data$Ileum_Dx == "positive"))
length(which((pcr.data$Flot %in% c("positive", "negative")) & !is.na(pcr.data$Ap5)))
length(which((pcr.data$Ap5 %in% c("positive", "negative")) & !is.na(pcr.data$Flot)))
length(which(pcr.data$Ap5 == "positive"))
length(which(pcr.data$Flot == "positive"))
length(which(pcr.data$Cecum_Dx == "positive"))
length(which(pcr.data$Ileum_Dx == "positive"))
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
rawdata <- makeMiceTable("/home/victor/Dokumente/Data_important/")

rawdata <- 
filt.data <- rawdata[names(rawdata) %in% c("Mouse_ID", "Transect", "Sex", "Longitude", "Latitude", "Year", "HI")]

##Set working directory LOCAL
setwd("/home/victor/Dokumente/Sequences/Samples_2017")

# Load haplogroups
haplo <- read.csv("Haplogroups_14_17.csv")
haplo <- haplo[-c(1:3,5:10)]
colnames(haplo) <- c("Mouse_ID", "Haplogroup")
haplo$Mouse_ID <- gsub(pattern = " ", replacement = "", x = haplo$Mouse_ID)
length(which(haplo$Haplogroup == "A"))
length(which(haplo$Haplogroup == "B"))
length(which(haplo$Haplogroup == "C"))

# Replace extra spaces
filt.data$Mouse_ID <- gsub(pattern = " ", replacement = "", x = filt.data$Mouse_ID)
pcr.data$Mouse_ID <- gsub(pattern = " ", replacement = "", x = pcr.data$Mouse_ID)
filt.data$Mouse_ID <- gsub(pattern = "Sk3173", replacement = "SK_3173", x = filt.data$Mouse_ID)

# Merge Data bases
Total.DB <- merge(filt.data, haplo , by="Mouse_ID")

Haplo.map(Total.DB, margin = .3)


HI.map(Total.DB)

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

ggplot(na.omit(Total.DB), aes(x = Haplogroup, y = HI)) + 
  geom_violin(aes(fill = Haplogroup), alpha = 0.7) + 
  scale_fill_manual(values = c("#FF3300","#66CC00", "#FFCC00"))+ geom_jitter(width = 0.1) + 
  labs(title="Distribution of Eimeria haplogroups according to host genotype",x="Haplogroup", y = "HI") + theme_bw() + 
  theme(plot.title = element_text(size=24, face = "bold"), axis.text= element_text(size = 16), axis.title=element_text(size=18,face="bold"))

ggplot(na.omit(Total.DB), aes(x = Transect, y = HI)) + 
  geom_violin(aes(fill = Transect), alpha = 0.7) + 
  scale_fill_manual(values = c("#FF3300","#66CC00", "#FFCC00"))+ geom_jitter(width = 0.1) + 
  labs(title="Distribution of HI according to transect",x="Transect", y = "HI") + theme_bw()

##Files for host species
setwd("/home/victor/Dokumente/Sequences/Samples_2017")
setwd("/home/victor/Dokumente/Sequences/Samples_2017/ORF470")

# Load host
host <- read.csv("Hosts.csv")
host <- read.csv("ORF470/Host_ORF.csv")


#######################HAPLOTYPE NETWORK GENERAL CODE ##############################################################
##Haplotype network 

## Written on the 30th August 2017 by Alice Balard/Victor Hugo Jarquin

#let's cheat, I have to git that back after:

# Function creating haplotypes:

##With ORF (Without any reference sequence)
VicAlign <- "/home/victor/Dokumente/Sequences/Samples_2017/Haplotype_network/ORF470_Haplotype_230218.fasta"
VicAlign <- "/home/victor/Dokumente/Sequences/Samples_2017/ORF470/ORF470_190418_Haplo_2.fasta"
VicAlign <- "ORF470/Concatenated_COI_ORF_Haplo.fasta"
##With 18S
VicAlign2 <- "/home/victor/Dokumente/Sequences/Samples_2017/Haplotype_network/18S_Haplotype_230218.fasta"

##With COI
VicAlign3 <- "/home/victor/Dokumente/Sequences/Samples_2017/Haplotype_network/COI_Haplotype_230218.fasta"
VicAlignM <- "Haplotype/COI_Rodents_Haplo_250518.fasta"

##COI_rodents
VicAlign8 <- "/home/victor/Dokumente/Sequences/Samples_2017/COI/COI_Rodent_Haplo_020318.fasta"

##Concatenated Haplotype
VicAlign4 <- "/home/victor/Dokumente/Sequences/Samples_2017/Haplotype_network/Concatenated_Haplotype_230218.fasta"
VicAlign5 <- "/home/victor/Dokumente/Sequences/Samples_2017/Haplotype_network/Concatenated_Haplotype_260218.fasta"
VicAlign6 <- "/home/victor/Dokumente/Sequences/Samples_2017/Haplotype_network/Concatenated_COI_ORF_270218.fasta"
VicAlign7 <- "/home/victor/Dokumente/Sequences/Samples_2017/Haplotype_network/Concatenated_18S_ORF_270218.fasta"


#Previous versions 

## ORF (E. falciformis/E. vermiformis/E. separata as references)
#VicAlign5 <- "/home/victor/Dokumente/Sequences/Samples_2014_2015/Haplotype_network/ORF470_Haplotype_310817.fasta"

## 18S (E. falciformis/E. vermiformis/E. ferrisi/E. separata as references)
#VicAlign7 <- "/home/victor/Dokumente/Sequences/Samples_2014_2015/Haplotype_network/n18S_Haplotype_310817.2.fasta"

## COI (E. falciformis/E. vermiformis/E. ferrisi/E. separata as references)
#VicAlign8 <- "/home/victor/Dokumente/Sequences/Samples_2014_2015/Haplotype_network/COI_Haplotype_310817.fasta"

##Concatenated Haplotype
#VicAlign5 <- "/home/victor/Dokumente/Sequences/Samples_2014_2015/Haplotype_network/Concatenated_haplotype_300817.fasta"
#VicAlign9 <- "/home/victor/Dokumente/Sequences/Samples_2014_2015/Haplotype_network/Concatenated_Haplotype_310817.fasta"
#VicAlign10 <- "/home/victor/Dokumente/Sequences/Samples_2014_2015/Haplotype_network/Concatenated_Haplotype_060917.fasta"

##Concatenated Madrid Apicowplexa
#VicAlign11 <- "/home/victor/Dokumente/Sequences/Samples_2014_2015/Haplotype_network/Concatenated_Haplotype_060917.2.fasta"


##ORF
d <- ape::read.dna(VicAlign, format = "fasta")

##18S  
d <-ape::read.dna(VicAlign2, format = "fasta")

##COI
d <- ape::read.dna(VicAlign3, format = "fasta")
d <- ape::read.dna(VicAlign8, format = "fasta")
##Concatenated
d <- ape::read.dna(VicAlign4, format = "fasta")
d <- ape::read.dna(VicAlign5, format = "fasta")
d <- ape::read.dna(VicAlign6, format = "fasta")
d <- ape::read.dna(VicAlign7, format = "fasta")

#d <- ape::read.dna(VicAlign5, format = "fasta")
#d <- ape::read.dna(VicAlign9, format = "fasta")
#d <- ape::read.dna(VicAlign10, format = "fasta")
# Apicoxplexa version d <- ape::read.dna(VicAlign11, format = "fasta")

## Turn the row names into corresponding HI :

seqnames <- data.frame(Mouse_ID = labels(d))
seqnames
matchname <- merge(seqnames, filt.data, by = "Mouse_ID", all.x = TRUE, sort = FALSE)
matchname
matchname$HI <- round(matchname$HI, 2)
matchname

## Turn the row names into corresponding Haplotype assignment 
seqnames <- data.frame(Mouse_ID = labels(d))
seqnames
matchname <- merge(seqnames, haplo, by = "Mouse_ID", all.x = TRUE, sort = FALSE)
matchname

## Turn the row names into corresponding host_species :

seqnames <- data.frame(Seq_ID = labels(d))
seqnames
matchname <- merge(seqnames, host, by = "Seq_ID", all.x = TRUE, sort = FALSE)
matchname

# Give ref names for ref sequences (Just in case Reference sequences are)
#matchname[grep("E_", matchname$Mouse_ID), ]$HI <- as.character(matchname[grep("E_", matchname$Mouse_ID), ]$Mouse_ID)

# Give name for NA data
matchname[which(is.na(matchname$HI)),]$HI <- "unknow_yet"
matchname
# Merge all
match <- merge(seqnames, matchname, by = "Mouse_ID", sort = FALSE)
match
## Change name of the sequence by corresponding HI :
rownames(d) <- match$HI

## Change name of the sequences by corresponding Haplogroup: 
rownames(d) <- match$Haplogroup

## Change name of the sequences by corresponding Host_species: 
rownames(d) <- matchname$Host_species

## Run haplotype analysis 

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

#mycols <- c(colorRampPalette(c("blue", "red"))(ncol(ind.hap) - 1), "grey") # "green", "darkgreen", "darkorange", "yellow", "hotpink",  "grey")


## Colors host_species 

mycols <- c("#91ED97","#03BA0F", "#015B07", "#B0038D", "#F02009", "#F78667", "#F86E06", "#1703F9", "#F8E1AB", "#F0E807", "#F781F3") # "green", "darkgreen", "darkorange", "yellow", "hotpink",  "grey")

# if you want to have the number of sequences in each circle:
attr(net, "labels") <- attr(net, "freq")

plot(net, size = .6*sqrt(attr(net, "freq")), pie = ind.hap, fast = T, scale.ratio = .75, 
     cex = 1, bg = mycols, lwd = 1, show.mutation = 2) 
legend(0, -10, colnames(ind.hap), fill = mycols, pch=NA, ncol=2, cex = .8) #+ scale_color_gradient("Hybrid\nindex", high="red",low="blue")


##ORF Haplo

##ORF Haplo colors
mycols <- c("#91ED97", "#015B07", "#CEA585", "#B0038D", "#1703F9","#F8573A", "#F0E807") 


plot(net, size = .7*sqrt(attr(net, "freq")), pie = ind.hap, scale.ratio = .7, 
     cex = 1, bg = mycols, show.mutation = 1) 
legend(-0, 25, colnames(ind.hap), fill = mycols, pch=NA, ncol=2, cex = 0.8)


##COI+ORF Haplo
mycols <- c("#91ED97", "#015B07", "#B0038D", "#1703F9", "#F0E807") 
attr(net, "labels") <- attr(net, "freq")

plot(net, size = .7*sqrt(attr(net, "freq")), pie = ind.hap, fast = TRUE, scale.ratio = 10, 
     cex = 1, bg = mycols, show.mutation = 1) 
legend(-0, -25, colnames(ind.hap), fill = mycols, pch=NA, ncol=2, cex = 0.8)


#+ scale_color_gradient("Hybrid\nindex", high="red",low="blue")
#}

# the size of the circles are proportinal to the number of sequences present in each haplotype (square root transformed)

## Colors haplotype 

#mycols <- c("#FF3300","#66CC00", "#FFCC00") # "green", "darkgreen", "darkorange", "yellow", "hotpink",  "grey")

#plot(net, size = attr(net, "freq"), pie = ind.hap, fast = TRUE, scale.ratio = 1, 
#     cex = 1, bg = mycols) 
#legend(55, -60, colnames(ind.hap), fill = mycols, pch=3, ncol=2, cex = 0.6) #+ scale_color_gradient("Hybrid\nindex", high="red",low="blue")
#}


## Colors host_species 

#mycols <- c("palegreen","limegreen", "seagreen", "violet", "salmon2", "red", "springgreen", "blue", "wheat", "yellow") # "green", "darkgreen", "darkorange", "yellow", "hotpink",  "grey")

#plot(net, size = attr(net, "freq"), pie = ind.hap, fast = TRUE, scale.ratio = 1, 
    # cex = 1, bg = mycols) 
#legend(-120, 20, colnames(ind.hap), fill = mycols, pch=3, ncol=2, cex = 0.6) #+ scale_color_gradient("Hybrid\nindex", high="red",low="blue")
#}

################# HAPLOTYPE NETWORK COI MANUSCRIPT #############

###Manuscript Haplotype
setwd("/home/victor/Dokumente/Sequences/Manuscript/")
host <- read.csv("Haplotype/Host_haplo.csv")
VicAlignM <- "Haplotype/COI_Rodents_Haplo_250518.fasta"
VicAlignM <- "Haplotype/COI_Rodents_Haplo_090718.fasta" ##Including Tabea's sequences
VicAlignM <- "Haplotype/COI_Rodents_Haplo_150818.fasta" ##Including last COI sequences
d <- ape::read.dna(VicAlignM, format = "fasta")
seqnames <- data.frame(Seq_ID = labels(d))
seqnames
matchname <- merge(seqnames, host, by = "Seq_ID", all.x = TRUE, sort = FALSE)
matchname
rownames(d) <- matchname$Host_species
rownames(d) <- matchname$Host_genus #Tabea#s preliminary host identification
e <- dist.dna(d)

h <- pegas::haplotype(d)

h <- sort(h, what = "label")

(net <- pegas::haploNet(h))

ind.hap<-with(
  stack(setNames(attr(h, "index"), rownames(h))),
  table(hap=ind, pop=rownames(d)[values])
)

mycols <- c("#91ED97","#03BA0F", "#015B07", "#B0038D", "#F02009", "#F78667", "#F86E06", "#1703F9", "#F8E1AB", "#F0E807", "#F781F3") # "green", "darkgreen", "darkorange", "yellow", "hotpink",  "grey")
mycols <- c("#03BA0F", "#F02009", "#1703F9", "#F8E1AB") #Tabeas haplotype
attr(net, "labels") <- attr(net, "freq")

plot(net, size = .6*sqrt(attr(net, "freq")), pie = ind.hap, fast = T, scale.ratio = .75, 
     cex = 1, bg = mycols, lwd = 1, show.mutation = 2) 
legend(0, -10, colnames(ind.hap), fill = mycols, pch=NA, ncol=2, cex = .8) #+ scale_color_gradient("Hybrid\nindex", high="red",low="blue")

########### BARCODE GAP ANALYSIS #############


COI <- ape::read.dna(VicAlign8, format = "fasta")


barcoding.gap(COI, dist="raw")
barcoding.gap(COI, dist="euclidean")
barcoding.gap(COI, dist="K80")

FMFtheta12(COI)

##Not finished

#Calculation intraspecific variation (sd) of the potential species theta1, and mean interspecific distance
#(here, the mean distance between the potential species and its nearest neighbor theta2) 
#(fuzzy-set based method,slightly modified from Zhang et al. 2012). The calculation was done for all species
#in the reference dataset.


##########STATISTICS GLM MANUSCRIPT########### 
###Regresion Model 
##Aim: Explain how likely is to amplify one of the 3 markers due to the fact of being assigned to an especific Eimeria species
#Analysis just takking into account COI with cocci primers
setwd("/home/victor/Dokumente/Sequences/Manuscript/")
phylogroup<- read.csv("Genotyped_mice.csv")
phylogroup$Mouse_ID <- gsub(pattern = " ", replacement = "", x = phylogroup$Mouse_ID)
phylogroup$Species <- gsub(pattern = " ", replacement = "", x = phylogroup$Species)
phylogroup$E_ferrisi <- gsub(pattern = " ", replacement = "", x = phylogroup$E_ferrisi)
phylogroup$E_vermiformis <- gsub(pattern = " ", replacement = "", x = phylogroup$E_vermiformis)

phylogroup<- phylogroup[phylogroup$Species%in%c("E_falciformis", "E_ferrisi", "E_vermiformis"),]
phylogroup$Species <- as.factor(as.character(phylogroup$Species))
table(phylogroup$Species) ##Single Infections in genotyped mice


##Probability to get a sequence from X genetic marker due to the fact to be assigned to a certain Eimeria species
##Considering all sequences from single infections N= 153 

logReg18S <- glm(n18S_Seq%in%"positive"~Species, family="binomial", data=phylogroup)
summary(logReg18S)

table(phylogroup$n18S_Seq%in%"positive", phylogroup$Species)

logRegCOI <- glm(COI_Seq%in%"positive"~Species, family="binomial", data=phylogroup)
summary(logRegCOI)

table(phylogroup$COI_Seq%in%"positive", phylogroup$Species)

logRegORF <- glm(ORF470_Seq%in%"positive"~Species, family="binomial", data=phylogroup)
summary(logRegORF)

table(phylogroup$ORF470_Seq%in%"positive", phylogroup$Species)


##posthoc 
##Multicomparison

##18S
glht(logReg18S, linfct = mcp(Species = "Tukey"))
summary(glht(logReg18S, linfct = mcp(Species = "Tukey")))

##COI
glht(logRegCOI, linfct = mcp(Species = "Tukey"))
summary(glht(logRegCOI, linfct = mcp(Species = "Tukey")))

##ORF
glht(logRegORF, linfct = mcp(Species = "Tukey"))
summary(glht(logRegORF, linfct = mcp(Species = "Tukey")))

##Graphs (Venn diagram by species)
# Venn diagram
Efalciformis <- phylogroup[phylogroup$Species%in%"E_falciformis",]
Eferrisi <- phylogroup[phylogroup$Species%in%"E_ferrisi",]
Evermiformis <- phylogroup[phylogroup$Species%in%"E_vermiformis",]
  
grid.newpage()
draw.triple.venn(area1 = nrow(subset(Efalciformis, n18S_Seq == "positive")), area2 = nrow(subset(Efalciformis, COI_Seq == "positive")), 
                 area3 = nrow(subset(Efalciformis, ORF470_Seq == "positive")), n12 = nrow(subset(Efalciformis, n18S_Seq == "positive"&COI_Seq== "positive")), 
                 n23 = nrow(subset(Efalciformis, COI_Seq== "positive"&ORF470_Seq =="positive")), n13 = nrow(subset(Efalciformis, n18S_Seq == "positive"&ORF470_Seq == "positive")), 
                 n123 = nrow(subset(Efalciformis, n18S_Seq == "positive"&ORF470_Seq == "positive"&COI_Seq=="positive")), category = c("18S", "COI", "ORF470"), 
                 lty = rep(1,3), col = c("dodgerblue4", "firebrick3", "darkgreen"), lwd = rep(2,3),
                 fill = c("dodgerblue4", "firebrick3", "darkgreen"), alpha = c(0.3, 0.3, 0.3), cex =2, cat.cex = 2.5, cat.default.pos = 'outer', 
                 cat.col = c("dodgerblue4", "firebrick3", "darkgreen"))

grid.newpage()
draw.triple.venn(area1 = nrow(subset(Eferrisi, n18S_Seq == "positive")), area2 = nrow(subset(Eferrisi, COI_Seq == "positive")), 
                 area3 = nrow(subset(Eferrisi, ORF470_Seq == "positive")), n12 = nrow(subset(Eferrisi, n18S_Seq == "positive"&COI_Seq== "positive")), 
                 n23 = nrow(subset(Eferrisi, COI_Seq== "positive"&ORF470_Seq =="positive")), n13 = nrow(subset(Eferrisi, n18S_Seq == "positive"&ORF470_Seq == "positive")), 
                 n123 = nrow(subset(Eferrisi, n18S_Seq == "positive"&ORF470_Seq == "positive"&COI_Seq=="positive")), category = c("18S", "COI", "ORF470"), 
                 lty = rep(1,3), col = c("dodgerblue4", "firebrick3", "darkgreen"), lwd = rep(2,3),
                 fill = c("dodgerblue4", "firebrick3", "darkgreen"), alpha = c(0.3, 0.3, 0.3), cex =2, cat.cex = 2.5, cat.default.pos = 'outer', 
                 cat.col = c("dodgerblue4", "firebrick3", "darkgreen"))

grid.newpage()
draw.triple.venn(area1 = nrow(subset(Evermiformis, n18S_Seq == "positive")), area2 = nrow(subset(Evermiformis, COI_Seq == "positive")), 
                 area3 = nrow(subset(Evermiformis, ORF470_Seq == "positive")), n12 = nrow(subset(Evermiformis, n18S_Seq == "positive"&COI_Seq== "positive")), 
                 n23 = nrow(subset(Evermiformis, COI_Seq== "positive"&ORF470_Seq =="positive")), n13 = nrow(subset(Evermiformis, n18S_Seq == "positive"&ORF470_Seq == "positive")), 
                 n123 = nrow(subset(Evermiformis, n18S_Seq == "positive"&ORF470_Seq == "positive"&COI_Seq=="positive")), category = c("18S", "COI", "ORF470"), 
                 lty = rep(1,3), col = c("dodgerblue4", "firebrick3", "darkgreen"), lwd = rep(2,3),
                 fill = c("dodgerblue4", "firebrick3", "darkgreen"), alpha = c(0.3, 0.3, 0.3), cex =2, cat.cex = 2.5, cat.default.pos = 'outer', 
                 cat.col = c("dodgerblue4", "firebrick3", "darkgreen"))

############Analysis considering COI get using cocci + Eim primers#######

phylogroup2<- read.csv("Genotyped_mice2.csv")
phylogroup2$Mouse_ID <- gsub(pattern = " ", replacement = "", x = phylogroup2$Mouse_ID)
phylogroup2$Species <- gsub(pattern = " ", replacement = "", x = phylogroup2$Species)
phylogroup2$E_ferrisi <- gsub(pattern = " ", replacement = "", x = phylogroup2$E_ferrisi)
phylogroup2$E_vermiformis <- gsub(pattern = " ", replacement = "", x = phylogroup2$E_vermiformis)

phylogroup2<- phylogroup2[phylogroup2$Species%in%c("E_falciformis", "E_ferrisi", "E_vermiformis"),]
phylogroup2$Species <- as.factor(as.character(phylogroup2$Species))
table(phylogroup2$Species) ##Single Infections in genotyped mice

###Data exploration by marker combination
length(which(phylogroup2$n18S_Seq == "positive" & phylogroup2$COI_Seq == "positive"))
length(which(phylogroup2$ORF470_Seq == "positive" & phylogroup2$COI_Seq == "positive"))
length(which(phylogroup2$n18S_Seq == "positive" & phylogroup2$ORF470_Seq == "positive"))

###By Species
length(which(phylogroup2$n18S_Seq == "positive" & phylogroup2$COI_Seq == "positive" & phylogroup2$E_falciformis == "positive"))
length(which(phylogroup2$n18S_Seq == "positive" & phylogroup2$COI_Seq == "positive" & phylogroup2$E_ferrisi == "positive"))
length(which(phylogroup2$n18S_Seq == "positive" & phylogroup2$COI_Seq == "positive" & phylogroup2$E_vermiformis == "positive"))

length(which(phylogroup2$n18S_Seq == "positive" & phylogroup2$ORF470_Seq == "positive" & phylogroup2$E_falciformis == "positive"))
length(which(phylogroup2$n18S_Seq == "positive" & phylogroup2$ORF470_Seq == "positive" & phylogroup2$E_ferrisi == "positive"))
length(which(phylogroup2$n18S_Seq == "positive" & phylogroup2$ORF470_Seq == "positive" & phylogroup2$E_vermiformis == "positive"))

length(which(phylogroup2$ORF470_Seq == "positive" & phylogroup2$COI_Seq == "positive" & phylogroup2$E_falciformis == "positive"))
length(which(phylogroup2$ORF470_Seq == "positive" & phylogroup2$COI_Seq == "positive" & phylogroup2$E_ferrisi == "positive"))
length(which(phylogroup2$ORF470_Seq == "positive" & phylogroup2$COI_Seq == "positive" & phylogroup2$E_vermiformis == "positive"))

logRegEfal <- glm(E_falciformis%in%"positive"~n18S_Seq+COI_Seq, family = "binomial", data = phylogroup2)
summary(logRegEfal)
table(phylogroup2$E_falciformis%in%"positive", phylogroup2$n18S_Seq, phylogroup2$COI_Seq)

##Probability to get a sequence from X genetic marker due to the fact to be assigned to a certain Eimeria species
##Considering all sequences from single infections N= 153 

logReg18S2 <- glm(n18S_Seq%in%"positive"~Species, family="binomial", data=phylogroup2)
summary(logReg18S2)
confint(logReg18S2)

table(phylogroup2$n18S_Seq%in%"positive", phylogroup2$Species)

logRegCOI2 <- glm(COI_Seq%in%"positive"~Species, family="binomial", data=phylogroup2)
summary(logRegCOI2)

table(phylogroup2$COI_Seq%in%"positive", phylogroup2$Species)

logRegORF2 <- glm(ORF470_Seq%in%"positive"~Species, family="binomial", data=phylogroup2)
summary(logRegORF2)

table(phylogroup2$ORF470_Seq%in%"positive", phylogroup2$Species)


##posthoc 
##Multicomparison

##18S
glht(logReg18S2, linfct = mcp(Species = "Tukey"))
summary(glht(logReg18S2, linfct = mcp(Species = "Tukey")))

##COI
glht(logRegCOI2, linfct = mcp(Species = "Tukey"))
summary(glht(logRegCOI2, linfct = mcp(Species = "Tukey")))

##ORF
glht(logRegORF2, linfct = mcp(Species = "Tukey"))
summary(glht(logRegORF2, linfct = mcp(Species = "Tukey")))

##Graphs (Venn diagram by species)
# Venn diagram
Efalciformis2 <- phylogroup2[phylogroup2$Species%in%"E_falciformis",]
Eferrisi2 <- phylogroup2[phylogroup2$Species%in%"E_ferrisi",]
Evermiformis2 <- phylogroup2[phylogroup2$Species%in%"E_vermiformis",]

grid.newpage()
draw.triple.venn(area1 = nrow(subset(Efalciformis2, n18S_Seq == "positive")), area2 = nrow(subset(Efalciformis2, COI_Seq == "positive")), 
                 area3 = nrow(subset(Efalciformis2, ORF470_Seq == "positive")), n12 = nrow(subset(Efalciformis2, n18S_Seq == "positive"&COI_Seq== "positive")), 
                 n23 = nrow(subset(Efalciformis2, COI_Seq== "positive"&ORF470_Seq =="positive")), n13 = nrow(subset(Efalciformis2, n18S_Seq == "positive"&ORF470_Seq == "positive")), 
                 n123 = nrow(subset(Efalciformis2, n18S_Seq == "positive"&ORF470_Seq == "positive"&COI_Seq=="positive")), category = c("18S", "COI", "ORF470"), 
                 lty = rep(1,3), col = c("dodgerblue4", "firebrick3", "darkgreen"), lwd = rep(2,3),
                 fill = c("dodgerblue4", "firebrick3", "darkgreen"), alpha = c(0.3, 0.3, 0.3), cex =2, cat.cex = 2.5, cat.default.pos = 'outer', 
                 cat.col = c("dodgerblue4", "firebrick3", "darkgreen"))

grid.newpage()
draw.triple.venn(area1 = nrow(subset(Eferrisi2, n18S_Seq == "positive")), area2 = nrow(subset(Eferrisi2, COI_Seq == "positive")), 
                 area3 = nrow(subset(Eferrisi2, ORF470_Seq == "positive")), n12 = nrow(subset(Eferrisi2, n18S_Seq == "positive"&COI_Seq== "positive")), 
                 n23 = nrow(subset(Eferrisi2, COI_Seq== "positive"&ORF470_Seq =="positive")), n13 = nrow(subset(Eferrisi2, n18S_Seq == "positive"&ORF470_Seq == "positive")), 
                 n123 = nrow(subset(Eferrisi2, n18S_Seq == "positive"&ORF470_Seq == "positive"&COI_Seq=="positive")), category = c("18S", "COI", "ORF470"), 
                 lty = rep(1,3), col = c("dodgerblue4", "firebrick3", "darkgreen"), lwd = rep(2,3),
                 fill = c("dodgerblue4", "firebrick3", "darkgreen"), alpha = c(0.3, 0.3, 0.3), cex =2, cat.cex = 2.5, cat.default.pos = 'outer', 
                 cat.col = c("dodgerblue4", "firebrick3", "darkgreen"))

grid.newpage()
draw.triple.venn(area1 = nrow(subset(Evermiformis2, n18S_Seq == "positive")), area2 = nrow(subset(Evermiformis2, COI_Seq == "positive")), 
                 area3 = nrow(subset(Evermiformis2, ORF470_Seq == "positive")), n12 = nrow(subset(Evermiformis2, n18S_Seq == "positive"&COI_Seq== "positive")), 
                 n23 = nrow(subset(Evermiformis2, COI_Seq== "positive"&ORF470_Seq =="positive")), n13 = nrow(subset(Evermiformis2, n18S_Seq == "positive"&ORF470_Seq == "positive")), 
                 n123 = nrow(subset(Evermiformis2, n18S_Seq == "positive"&ORF470_Seq == "positive"&COI_Seq=="positive")), category = c("18S", "COI", "ORF470"), 
                 lty = rep(1,3), col = c("dodgerblue4", "firebrick3", "darkgreen"), lwd = rep(2,3),
                 fill = c("dodgerblue4", "firebrick3", "darkgreen"), alpha = c(0.3, 0.3, 0.3), cex =2, cat.cex = 2.5, cat.default.pos = 'outer', 
                 cat.col = c("dodgerblue4", "firebrick3", "darkgreen"))



#############Association according to host #################

# Stacked Bar Plot with Colors and Legend
counts <- table(phylogroup$Species, phylogroup$Class)

barplot(counts, ylim = c(0,120), main="Eimeria species distribution by host type",
        xlab= "Host type", ylab= "Frequency", col=c("#FF3300","#66CC00", "#FFCC00"),
        legend = rownames(counts)) 


x1  = factor(counts, levels=c("Mmd", "Hyb", "Mmm"))
Counts.DB <- as.data.frame(counts) 
colnames(Counts.DB) <- c("Species", "Host", "Freq") 
Counts.DB

##Probability to be X Eimeria species due to the fact to be Y host type 
##Host assigned to 3 categories (just with genotyped mice N=139)

#E_falciformis
table(phylogroup$E_falciformis%in%"positive", phylogroup$Class)
logRegEfal <- glm(E_falciformis%in%"positive"~Class, family="binomial", data=phylogroup)
summary(logRegEfal)

#E_ferrisi
table(phylogroup$E_ferrisi%in%"positive", phylogroup$Class)
logRegEfer <- glm(E_ferrisi%in%"positive"~Class, family="binomial", data=phylogroup)
summary(logRegEfer)

#E_vermiformis
table(phylogroup$E_vermiformis%in%"positive", phylogroup$Class)
logRegEver <- glm(E_vermiformis%in%"positive"~Class, family="binomial", data=phylogroup)
summary(logRegEver)

##posthoc 
##Multicomparison

##E_falciformis
glht(logRegEfal, linfct = mcp(Class = "Tukey"))
summary(glht(logRegEfal, linfct = mcp(Class = "Tukey")))

##E_ferrisi
glht(logRegEfer, linfct = mcp(Class = "Tukey"))
summary(glht(logRegEfer, linfct = mcp(Class = "Tukey")))

##E_vermiformis
glht(logRegEver, linfct = mcp(Class = "Tukey"))
summary(glht(logRegEver, linfct = mcp(Class = "Tukey")))

##Plot GLM (check later)
#ggplot(phylogroup, aes(x=phylogroup$Class , y= E_falciformis)) + geom_point() + 
 # stat_smooth(method="glm", family="binomial", se=FALSE)

##Host genotype as continous variable 

summary(glm(E_falciformis%in%"positive" ~ HI, family="binomial", data=phylogroup)) # + violin plot -> more meaningful


##Violin plot for single infected mice
ggplot(na.omit(phylogroup), aes(x = Species, y = HI)) + 
  geom_violin(aes(fill = Species), alpha = 0.7) + 
  scale_fill_manual(values = c("#FF3300","#66CC00", "#FFCC00"))+ geom_jitter(width = 0.1) + 
  labs(title="Distribution of Eimeria species according to host genotype",x="Eimeria species", y = "Hybrid Index (HI)") + theme_bw() + 
  theme(plot.title = element_text(size=18, face = "bold"), axis.text= element_text(size = 16), axis.title=element_text(size=18,face="bold"))


#Table with HI and Host assignment


Total.res <- merge(filt.data, pcr.data, by="Mouse_ID")
length(which(Total.res$Ap5 %in% "negative"))


Total.DB <- read.csv("Total_DB.csv")
Total.DB$Mouse_ID <- gsub(pattern = " ", replacement = "", x = Total.DB$Mouse_ID)
Total.DB$Ap5 <- gsub(pattern = " ", replacement = "", x = Total.DB$Ap5)

finalData <- subset(x = Total.DB,
                    subset = (Total.DB$Ap5 %in% c("positive", "negative")))

length(which(finalData$Ap5 %in% "negative"))

######### OOCYST MEASURMENTS MANUSCRIPT#############

# Plot for compare measurements with references 
setwd("/home/victor/Dokumente/Oocysts/")
read.csv(file = "Comparative_measurments.csv")
measurments <- read.csv(file = "Comparative_measurments.csv")
plotConfidence(measurments$mean, lower =  measurments$lower, upper =  measurments$higer, xlim = c(0.9,1.6),
               labels = c("Eimeria falciformis wild", "Eimeria falciformis ref", "Eimeria ferrisi wild", "Eimeria falciformis ref","Eimeria vermiformis wild", "Eimeria falciformis ref"),
               lwd = 2, title.labels = "Strain", col = "black", refline = FALSE, rightmargin = 0.0001, leftmargin = 0.0001)

##Import data 
setwd("/home/victor/Dokumente/Oocysts/")
raw.data <- read.csv("Measurements_oocysts.csv")

oocyst.data <- raw.data[-c(2,9:12)]


names(oocyst.data)[names(oocyst.data) == "S_Lenght"] <- "S_Length"

####Plot comparison among isolates E. ferrisi
names(oocyst.data)[names(oocyst.data) == "Strain_ID"] <- "Mouse_ID"

levels(oocyst.data$Close.related)

oocyst.data <- merge(mouse.table, oocyst.data, all.y=TRUE)

ggplot(data= oocyst.data[oocyst.data$Close.related=="E. ferrisi",], aes(HI, S_LW_Ratio,Close.related)) + 
  geom_jitter() 

ggplot(oocyst.data[oocyst.data$Close.related=="E. ferrisi",],
       aes(x = Mouse_ID, y = O_LW_Ratio, col = HI)) +
  scale_color_gradient(low = "blue", high = "red") +
  geom_boxplot() +
  geom_point()

oocyst.data$Close.related

##Oocyst and sporocysts L/W ratio statistics
LWratioDF <- group_by(oocyst.data,Close.related) %>%
  dplyr::summarise(
    count = n(),
    meanO_LW_Ratio = mean(O_LW_Ratio, na.rm = TRUE),
    sdO_LW_Ratio = sd(O_LW_Ratio, na.rm = TRUE),
    seO_LW_Ratio = sdO_LW_Ratio / sqrt(count),  # Calculate standard error of the mean
    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult = qt(.95/2 + .5, n()-1),
    ciO_LW_Ratio = seO_LW_Ratio * ciMult,
    ciminO_LW_Ratio = meanO_LW_Ratio - ciO_LW_Ratio,
    cimaxO_LW_Ratio = meanO_LW_Ratio + ciO_LW_Ratio,
    ## and for Sporocysts
    meanS_LW_Ratio = mean(S_LW_Ratio, na.rm = TRUE),
    sdS_LW_Ratio = sd(S_LW_Ratio, na.rm = TRUE),
    seS_LW_Ratio = sdS_LW_Ratio / sqrt(count),  # Calculate standard error of the mean
    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult = qt(.95/2 + .5, n()-1),
    ciS_LW_Ratio = seS_LW_Ratio * ciMult,
    ciminS_LW_Ratio = meanS_LW_Ratio - ciS_LW_Ratio,
    cimaxS_LW_Ratio = meanS_LW_Ratio + ciS_LW_Ratio
  )

##one way ANOVA 
oocyst.anova <- aov(O_LW_Ratio ~ Close.related, data = oocyst.data)
summary(oocyst.anova)

Socyst.anova <- aov(S_LW_Ratio ~ Close.related, data = oocyst.data)
summary(Socyst.anova)

##posthoc 
TukeyHSD(oocyst.anova)
TukeyHSD(Socyst.anova)
##Multicomparison
glht(oocyst.anova, linfct = mcp(Close.related = "Tukey"))
summary(glht(oocyst.anova, linfct = mcp(Close.related = "Tukey")))


glht(Socyst.anova, linfct = mcp(Close.related = "Tukey"))
summary(glht(Socyst.anova, linfct = mcp(Close.related = "Tukey")))

##95% CI for measurments 
group.CI(O_LW_Ratio ~ Close.related, data = oocyst.data, ci = 0.95)

group.CI(S_LW_Ratio ~ Close.related, data = oocyst.data, ci = 0.95)


##Boxplot
fill <- c("#FF3300","#66CC00", "#FFCC00")
line <- "black"

# Choice plot type
victor_plot <- geom_boxplot(fill = fill, colour = line,
                            alpha = 0.8, outlier.colour = "black", outlier.shape = 20)

geom_boxplot

## L/W ratio Oocysts
plot_oocyst <- ggplot(oocyst.data, aes(x = Close.related, y = O_LW_Ratio)) +
  scale_x_discrete(name = "Eimeria species") +
  scale_y_continuous(name = "Oocyst length/width ratio",
                     limits=c(0.95, 1.60)) + ggtitle("Oocysts L/W ratio by group")+
  victor_plot + theme_bw() +
  theme(plot.title = element_text(size=24, face = "bold"), axis.text= element_text(size = 18), axis.title=element_text(size=18,face="bold")) +
  geom_jitter(pch= 21)


##L/W ratio Sporocysts 
plot_sporo <- ggplot(oocyst.data, aes(x = Close.related, y = S_LW_Ratio)) +
  scale_x_discrete(name = "Eimeria species") +
  scale_y_continuous(name = "Sporocyst length/width ratio",
                     limits=c(1.1, 2.35)) + ggtitle("Sporocysts L/W ratio by group")+
  victor_plot + theme_bw() +
  theme(plot.title = element_text(size=24, face = "bold"), axis.text= element_text(size = 18), axis.title=element_text(size=18,face="bold")) +
  geom_jitter(pch= 21)


grid.arrange(plot_oocyst, plot_sporo, ncol = 2)


#Mean plots Alice method with ggplot

oocyst.data <- merge(oocyst.data, LWratioDF)

##Oocyst 
ciOocyst <- ggplot(oocyst.data, aes(oocyst.data$Close.related, y = oocyst.data$O_LW_Ratio)) +
  geom_violin(aes(col = Close.related), trim = F, alpha = 0.8) +
  geom_jitter(aes(fill = Close.related),
              width=0.1, size=3, alpha = 0.4, pch = 21) +
  geom_point(aes(y = meanO_LW_Ratio), size = 2) +
  scale_color_manual(values = c("#FF3300","#66CC00", "#FFCC00")) +
  scale_fill_manual(values = c("#FF3300","#66CC00", "#FFCC00")) +
  geom_errorbar(aes(ymin = ciminO_LW_Ratio, ymax = cimaxO_LW_Ratio), size = 0.80, width=0.25) +
  labs(y = "Oocyst length/width ratio", x = "Eimeria species") + #ggtitle("Oocysts L/W ratio by group") +
  theme(legend.position="none", axis.title = element_text(size=18,face="bold"), axis.text = element_text(size = 18)) +
  theme_bw()


#Sporocyst
ciSporo <- ggplot(oocyst.data, aes(oocyst.data$Close.related, y = oocyst.data$S_LW_Ratio)) +
  geom_violin(aes(col = Close.related), trim = F, alpha = 0.8) +
  geom_jitter(aes(fill = Close.related),
              width=0.1, size=3, alpha = 0.4, pch = 21) +
  guides(fill=FALSE) + scale_fill_discrete(guide=FALSE) +
  geom_point(aes(y = meanS_LW_Ratio), size = 2) +
  scale_color_manual(values = c("#FF3300","#66CC00", "#FFCC00")) +
  scale_fill_manual(values = c("#FF3300","#66CC00", "#FFCC00")) +
  geom_errorbar(aes(ymin = ciminS_LW_Ratio, ymax = cimaxS_LW_Ratio), size = 0.80, width=0.25) +
  labs(y = "Sporocyst length/width ratio", x = "Eimeria species") + #ggtitle("Sporocysts L/W ratio by group") +
  theme(legend.position="none", axis.title = element_text(size=18,face="bold"), axis.text = element_text(size = 18)) +
  theme_bw()


grid.arrange(ciOocyst, ciSporo, ncol = 2)

##Other options

# Mean plots
# ++++++++++++++++++++
# Plot weight by group
# Add error bars: mean_se
# Add confidence intervals: mean_ci
# Add Standar deviation: mean_sd
# (other values include: median_iqr, ....)

ggline(oocyst.data, x = "Close.related", y = "O_LW_Ratio", 
       add = c("mean_ci", "violin", "jitter"), point.color = "Black",
       ylab = "Oocyst length/width ratio", xlab = "Eimeria species", 
       color = "Close.related", palette = c("#FF3300","#66CC00", "#FFCC00"))

# Box plots
# ++++++++++++++++++++
# Plot weight by group and color by group

ggboxplot(oocyst.data, x = "Close.related", y = "O_LW_Ratio", 
          color = "Close.related", palette = c("#FF3300","#66CC00", "#FFCC00"),
          ylab = "Oocyst length/width ratio", xlab = "Eimeria species")


######qPCR 2016#######

rawqpcr <- read.csv("/home/victor/Dokumente/qPCR/Mouse_2016/qPCR_2016.csv")

ILvsCe <- ggplot(rawqpcr, aes(x= rawqpcr$CE_Neg_Dct, y= rawqpcr$IL_Neg_Dct, color=rawqpcr$Eimeria_ID)) +
  geom_jitter(shape=16, position=position_jitter(0.2)) + geom_hline(yintercept=-6, linetype="dashed", color = "red") + 
  geom_vline(xintercept = -6, linetype="dashed", color = "blue") + 
  labs(title="Effect of ileum infection on cecum detection",x="- ΔCt Cecum", y = "- ΔCt Ileum") +
 # geom_density_2d() +
  theme_bw()+
  theme(plot.title = element_text(color = "black", size = 22, face = "bold"), axis.title.x = element_text(color = "black", size = 18, face = "bold"), 
        axis.title.y = element_text(color = "black", size = 18, face = "bold"),
        axis.text.x = element_text(color = "dark grey", size = 16),
        axis.text.y = element_text(color = "dark grey", size = 16), 
        axis.line = element_line(size = 1), legend.text = element_text(size = 14))
  #theme_classic()

CevsIL <- ggplot(rawqpcr, aes(x= rawqpcr$IL_Neg_Dct, y= rawqpcr$CE_Neg_Dct,  color=rawqpcr$Eimeria_ID)) +
  geom_jitter(shape=16, position=position_jitter(0.2)) + geom_hline(yintercept=-6, linetype="dashed", color = "red") + 
  geom_vline(xintercept = -6, linetype="dashed", color = "blue") + 
  labs(title="Effect of cecum infection on ileum detection",x="- ΔCt Ileum", y = "- ΔCt Cecum") +
  theme_classic()

###Lets merge rawqpcr with eimeria species

rawEimeria <- read.csv("/home/victor/Dokumente/Sequences/Manuscript/Species_assignment_14_17.csv")

rawEimeria$Mouse_ID <- gsub(pattern = " ", replacement = "", x = rawEimeria$Mouse_ID)

species <- rawEimeria[names(rawEimeria) %in% c("Mouse_ID", "Assignment")]

qpcr.species <- merge(rawqpcr, species , by="Mouse_ID")


##Plot with Eimeria species 

ggplot(qpcr.species, aes(x= qpcr.species$CE_Neg_Dct, y= qpcr.species$IL_Neg_Dct, color=qpcr.species$Assignment, shape=qpcr.species$Eimeria_ID)) + 
  #geom_jitter(position=position_jitter(0.2)) + 
  geom_point()+
  #geom_text(label=rownames(qpcr.species$Mouse_ID)) +
  geom_hline(yintercept=-6, linetype="dashed", color = "black") + 
  geom_vline(xintercept = -6, linetype="dashed", color = "black") + 
  labs(title="Eimeria detection in ileum and cecum",x="- ΔCt Cecum", y = "- ΔCt Ileum") +
  scale_shape_manual(values=c(16, 17))+ 
  scale_color_manual(values=c("#FF00FF","#FF3300", "#66CC00", "#FFCC00", "#1703F9"))+
  scale_size_manual(values= 40)+
  #geom_density_2d() +
  theme_bw()+
  theme(plot.title = element_text(color = "black", size = 22, face = "bold"), axis.title.x = element_text(color = "black", size = 18, face = "bold"), 
        axis.title.y = element_text(color = "black", size = 18, face = "bold"),
        axis.text.x = element_text(color = "dark grey", size = 16),
        axis.text.y = element_text(color = "dark grey", size = 16), 
        axis.line = element_line(size = 1), legend.text = element_text(size = 14))


###With treshold at -5 XD
ggplot(qpcr.species, aes(x= qpcr.species$CE_Neg_Dct, y= qpcr.species$IL_Neg_Dct, color=qpcr.species$Assignment, shape=qpcr.species$Eimeria_ID)) + 
  #geom_jitter(position=position_jitter(0.2)) + 
  geom_point()+
  #geom_text(label=rownames(qpcr.species$Mouse_ID)) +
  geom_hline(yintercept=-5, linetype="dashed", color = "black") + 
  geom_vline(xintercept = -5, linetype="dashed", color = "black") + 
  labs(title="Eimeria detection in ileum and cecum",x="- ΔCt Cecum", y = "- ΔCt Ileum") +
  scale_shape_manual(values=c(16, 17))+ 
  scale_color_manual(values=c("#FF00FF","#FF3300", "#66CC00", "#FFCC00", "#1703F9"))+
  scale_size_manual(values= 40)+
  #geom_density_2d() +
  theme_bw()+
  theme(plot.title = element_text(color = "black", size = 22, face = "bold"), axis.title.x = element_text(color = "black", size = 18, face = "bold"), 
        axis.title.y = element_text(color = "black", size = 18, face = "bold"),
        axis.text.x = element_text(color = "dark grey", size = 16),
        axis.text.y = element_text(color = "dark grey", size = 16), 
        axis.line = element_line(size = 1), legend.text = element_text(size = 14))



        #scale_color_brewer(palette="Dark2")
#theme_classic()
#pdf(file = "/home/victor/Dokumente/Sequences/Manuscript/Figures/Tissue_qpcr.pdf", width = 9, height = 10, paper = "a4")
#1703F9 Negative
##FF3300 Eimeria falciformis 
#66CC00 Eimeria ferrisi
#FFCC00 Eimeria vermiformis
#FF00FF Double

##Venn diagram 

grid.newpage()
draw.pairwise.venn(area1 = nrow(subset(qpcr.species, qpcr.species$CE_Status == "Positive")), area2 = nrow(subset(qpcr.species, qpcr.species$IL_Status == "Positive")), cross.area = nrow(subset(qpcr.species, qpcr.species$CE_Status == "Positive"&qpcr.species$IL_Status == "Positive")), 
                   category = c("Cecum Dx", "Ileum Dx"), lty = rep("blank", 2), fill = c("light blue", "light green"), 
                   alpha = rep(0.5, 2), cat.pos = c(0, 180))


draw.triple.venn(area1 = nrow(subset(qpcr.species, n18S_Seq == "positive")), area2 = nrow(subset(tissue.data, COI_Seq == "positive")), 
                 area3 = nrow(subset(qpcr.species, ORF470_Seq == "positive")), n12 = nrow(subset(tissue.data, n18S_Seq == "positive"&COI_Seq== "positive")), 
                 n23 = nrow(subset(qpcr.species, COI_Seq== "positive"&ORF470_Seq =="positive")), n13 = nrow(subset(tissue.data, n18S_Seq == "positive"&ORF470_Seq == "positive")), 
                 n123 = nrow(subset(qpcr.species, n18S_Seq == "positive"&ORF470_Seq == "positive"&COI_Seq=="positive")), category = c("18S", "COI", "ORF470"), 
                 lty = rep(1,3), col = c("dodgerblue4", "firebrick3", "darkgreen"), lwd = rep(2,3),
                 fill = c("dodgerblue4", "firebrick3", "darkgreen"), alpha = c(0.3, 0.3, 0.3), cex =2, cat.cex = 2.5, cat.default.pos = 'outer', 
                 cat.col = c("dodgerblue4", "firebrick3", "darkgreen"))

