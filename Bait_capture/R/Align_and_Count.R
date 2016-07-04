###libraries
setwd("/SAN/Alices_sandpit/sequencing_data_dereplicated/")
#source("http://bioconductor.org/biocLite.R")
#biocLite("BiocUpgrade")
library(Rsubread)
library("rbamtools")
library(IRanges)
library(GenomicRanges)
library(Rsamtools)
library(ggplot2)
library(xtable)
library(compare)
library(ggrepel)
library(grid)
library(RColorBrewer)

themeAlice<-theme(text = element_text(size=15), 
                  #axis.line.x = element_line(linetype = "solid"),
                  axis.line.y = element_line(linetype = "solid"),
                  axis.line.x = element_line(linetype = "solid"),
                  axis.title = element_text(size=20),
                  axis.text = element_text(size=16),
                  axis.title.y=element_text(margin=margin(0,20,0,0)),
                  axis.title.x=element_text(margin=margin(20,0,0,0)),
                  panel.background = element_blank())
#legend.position = "none")


##### PART I





#Rsubread package: (1) build an index (E.Falciformis genome)

##################
##Index building##
##################
#An index needs to be built before read mapping can be performed. 
#buildindex creates a hash table for the reference genome, 
#which can then be used by Subread and Subjunc aligners for read alignment:
buildindex(basename="reference_index",reference="Efal_genome.fa")
#Rsubread creates a hash table for indexing the reference genome



#Tally program on Bash (2): Deduplication of sequence fragments [Tally processes both files record-by-record and pair up records at the same offset. This requires the option --pair-by-offset.]



#Rsubread package: (3) align the dereplicated sequences to the index
#Map paired-end reads
FilesR1 <- list.files(pattern = "R1_001.fastq.unique.gz")
FilesR2 <-list.files(pattern = "R2_001.fastq.unique.gz")
#secondList <- list(rep("sequencing_data",length(FilesR1)))
FilesR1path <- NA
FilesR2path <- NA

for (i in 1:length(FilesR1))
  (FilesR1path <- c(FilesR1path,paste(c(FilesR1[i]),collapse = "")))
FilesR1path <- FilesR1path[-1]
for (i in 1:length(FilesR2))
  (FilesR2path <- c(FilesR2path,paste(c(FilesR2[i]),collapse = "")))
FilesR2path <- FilesR2path[-1]
FilesR1path;FilesR2path

#Create 1 alignments (20min/alignment !!)
#*************
####for (i in 1:length(FilesR1path)){
####  reads1 <- FilesR1path[i]
####  reads2 <- FilesR2path[i]
####  Rsubread::align(index="reference_index",readfile1=reads1,readfile2=reads2, type="dna",maxMismatches=20, indels=10,
####                  output_file = paste("Alignment_", substr(FilesR1[i],1,6), sep = ""))
####}
##########
# repair manually the errors (to delete later) REPEATED FAILURE
####reads1 <- FilesR1path[9]
####reads2 <- FilesR2path[9]
####Rsubread::align(index="reference_index",readfile1=reads1,readfile2=reads2, type="dna",maxMismatches=20, indels=10,
####                  output_file = paste("Alignment_", substr(FilesR1[9],1,6), sep = ""))

####reads1 <- FilesR1path[11]
####reads2 <- FilesR2path[11]
####Rsubread::align(index="reference_index",readfile1=reads1,readfile2=reads2, type="dna",maxMismatches=20, indels=10,
####                  output_file = paste("Alignment_", substr(FilesR1[11],1,6), sep = ""))

####reads1 <- FilesR1path[15]
####reads2 <- FilesR2path[15]
####Rsubread::align(index="reference_index",readfile1=reads1,readfile2=reads2, type="dna",maxMismatches=20, indels=10,
####                  output_file = paste("Alignment_", substr(FilesR1[15],1,6), sep = ""))
#*************


######################
##Mapping percentage##
######################
#Function propmapped returns the proportion of mapped reads included in a SAM/BAM
#file. For paired end reads, it can return the proportion of mapped fragments (ie. read pairs).
#*************
for (i in c(1:8,10,12:14,16:19)) {
  assign(paste("PropMapped.dereplicated_", substr(FilesR1[i],1,6), sep = ""),
         propmapped(paste("Alignment_", substr(FilesR1[i],1,6), sep = "")))
}
#*************
# Add them all in ONE dataframe, make a distribution plot
#*************
name <- NA
for (i in c(1:8,10,12:14,16:19)) {
  name <- c(name,paste(substr(FilesR1[i],1,6)))
}
name <- name[-1]
#*************
Prop <- NA
for (i in c(1:8,10,12:14,16:19)) {
  prop <- paste("PropMapped.dereplicated_", substr(FilesR1[i],1,6), sep = "")
  prop1 <- get(prop) #get a list of values
  prop2 <- prop1$PropMapped #get the value
  Prop <- c(Prop, prop2)
}
Prop <- Prop[-1]

#*************
dat <- data.frame(name = name, propMapped = Prop)
dat$propMapped <- dat$propMapped*100
#*************

# Histogram overlaid with kernel density curve
Histo <- ggplot(dat, aes(x=propMapped)) + 
  geom_histogram(colour="black", fill="white")+ # Histogram with count on y-axis
  themeAlice+
  scale_x_continuous(name = "Proportion Mapped")+
  scale_y_continuous(breaks = 0:10)
Histo
#ggsave("PlotHistoDens.pdf")
#*************
x = knitr::kable(dat)
pdf("PropMapped.pdf", height=11, width=8.5)
grid.table(x)
dev.off()

###############################################
## Counting mapped reads for genomic features##
###############################################
fc2 <- list()
for (i in  c(1:8,10,12:14,16:19)) {
  fc2[[FilesR1[i]]] <-  featureCounts(paste("Alignment_", substr(FilesR1[i],1,6), sep = ""),
                                      annot.ext="MYbaits_Eimeria_V1.single120_feature_counted.gtf",
                                      isGTFAnnotationFile=TRUE,
                                      GTF.featureType = "sequence_feature",
                                      useMetaFeatures=FALSE,
                                      GTF.attrType="bait_id",
                                      # parameters specific to paired end reads
                                      isPairedEnd=TRUE,
                                      reportReads=TRUE)
}

mylist <- lapply(fc2, "[[", "counts")

mylist2 <- do.call(cbind,mylist)

head(mylist2)

table(rowSums(mylist2)>100)

table(rowSums(mylist2>10)>5) #a coverage >10 in >5 lib = how many baits are "working" with more than 10 sequences captured in more than 5 libraries

table(rowSums(mylist2>0)>10)#a coverage >0 in >10 lib = how many baits are "working" in more than 10 libraries



##### 

##### PART II
#PropDataFrame <- propmapped(paste("Alignment_", substr(FilesR1[1],1,6), sep = ""))
#for (i in c(2:8,10,12:14,16:19)) {
  PropDataFrame <- rbind(PropDataFrame, propmapped(paste("Alignment_", substr(FilesR1[i],1,6), sep = "")))
}


# Number of mapped reads/fragments will be counted and fraction of such reads/fragments will be calculated.


# Save as PDF
PropDataFrame
library(gridExtra)
pdf("PropDataFrame.pdf", height=11, width=8.5)
grid.table(PropDataFrame)
dev.off()
#

# :) Shorten the names for better visualisation
PropDataFrame$Samples<-base::gsub("Alignment_","",PropDataFrame$Samples)

newdata <- PropDataFrame[base::order(PropDataFrame$Samples),]
#add a column to see if single, double, digested
for (i in 1:nrow(newdata)){
  if (grepl("Do",newdata$Samples[i])) {
    newdata$type[i] <- "Double"
  }
  else if (grepl("Si",newdata$Samples[i])) {
    newdata$type[i] <- "Single"
  }
  else if (grepl("Di",newdata$Samples[i])) {
    newdata$type[i] <- "Digested"
  }
}

#!!! "undetermine" is "2808Di"
newdata$Samples[16] <- "2808Di"
newdata$type[16] <- "Digested"

# Nbr in million
newdata$NumTotalinmillion <-newdata$NumTotal/(1000000)

myplot<-ggplot(newdata, aes(x=PropMapped, y=NumTotalinmillion, label=Samples, col=type))+
  scale_color_manual(values=c(3,2,1), name="")+
  geom_point()+
  geom_label_repel(aes(label=Samples),label.size =1, size=5)+
  labs(x="Proportion of contigs mapped to the genome", y="Total number of contigs (in millions)")+
  themeAlice
#annotation_custom(grob=circleGrob(r=unit(1,"npc"), gp=gpar(fill="orange", alpha=0.2)), xmin=0.161, xmax=0.179, ymin=10.1, ymax=15.9)+
#annotate("text", x = rep(0.17,3), y = 12:14, size=8, label = c("Digested","Double","Single"),col=c(3,2,1))
myplot

###########################################
#add a column to see Eastern/Western clade#
###########################################
#2672 Eastern musculus
#2807 Eastern domesticus
#2808 Western domesticus
#2809 Western domesticus
#2811 Eastern domesticus
#2812 Eastern domesticus
#2848 Eastern hybrid
#2919 Eastern domesticus
#2TRAnn X apodemus
#95Anna X apodemus
#f9Anna X apodemus

# East/west classif: 
newdata$clade <- ifelse(grepl("2672",newdata$Samples),"Eastern",NA)
newdata$clade <- ifelse(grepl("2807",newdata$Samples),"Eastern",newdata$clade)
newdata$clade <- ifelse(grepl("2808",newdata$Samples),"Western",newdata$clade)
newdata$clade <- ifelse(grepl("2809",newdata$Samples),"Western",newdata$clade)
newdata$clade <- ifelse(grepl("2811",newdata$Samples),"Eastern",newdata$clade)
newdata$clade <- ifelse(grepl("2812",newdata$Samples),"Eastern",newdata$clade)
newdata$clade <- ifelse(grepl("2848",newdata$Samples),"Eastern",newdata$clade)
newdata$clade <- ifelse(grepl("2919",newdata$Samples),"Eastern",newdata$clade)
newdata$clade <- ifelse(grepl("2TRAnn",newdata$Samples),"X",newdata$clade)
newdata$clade <- ifelse(grepl("95Anna",newdata$Samples),"X",newdata$clade)
newdata$clade <- ifelse(grepl("f9Anna",newdata$Samples),"X",newdata$clade)

# Mus/dom/Apo classif: 
newdata$mice <- ifelse(grepl("2672",newdata$Samples),"musculus",NA)
newdata$mice <- ifelse(grepl("2807",newdata$Samples),"domesticus",newdata$mice)
newdata$mice <- ifelse(grepl("2808",newdata$Samples),"domesticus",newdata$mice)
newdata$mice <- ifelse(grepl("2809",newdata$Samples),"domesticus",newdata$mice)
newdata$mice <- ifelse(grepl("2811",newdata$Samples),"domesticus",newdata$mice)
newdata$mice <- ifelse(grepl("2812",newdata$Samples),"domesticus",newdata$mice)
newdata$mice <- ifelse(grepl("2848",newdata$Samples),"hybrid",newdata$mice)
newdata$mice <- ifelse(grepl("2919",newdata$Samples),"domesticus",newdata$mice)
newdata$mice <- ifelse(grepl("2TRAnn",newdata$Samples),"apodemus",newdata$mice)
newdata$mice <- ifelse(grepl("95Anna",newdata$Samples),"apodemus",newdata$mice)
newdata$mice <- ifelse(grepl("f9Anna",newdata$Samples),"apodemus",newdata$mice)


myplot3<-ggplot(newdata, aes(x=PropMapped, y=NumTotalinmillion, label=Samples, col=clade))+
  scale_color_manual(values=c("green","darkgreen","orange"), name="")+
  annotation_custom(grob=circleGrob(r=unit(1,"npc"), gp=gpar(col="red", alpha=1)), xmin=0.02, xmax=0.04, ymin=17.5)+
  geom_point(size=3)+
  geom_label_repel(aes(label=Samples),label.size =1, size=5)+
  labs(x="Proportion of contigs mapped to the genome", y="Total number of contigs (in millions)")+
  themeAlice+
  annotate("text", x = rep(0.16,3), y = 12:14, size=8, label = c("Closely related to Western clade","Eastern clade","Western clade"),col=c("green","darkgreen","orange"))
myplot3

###### ggplot function
#Create a custom color scale
myColors <- brewer.pal(12,"Paired")
names(myColors) <- levels(newdata$pairs)
colScale <- scale_colour_manual(name = "pairs",values = myColors)
######

#add a column to see the pairs
#newdata$pairs<- c(1,2,2,2,3,3,4,4,4,5,6,6,7,8,8,9,10,11,12)
#newdata$pairs <- as.factor(newdata$pairs)

#Plot
#myplot2<-ggplot(newdata, aes(x=PropMapped, y=NumTotalinmillion, label=Samples, col=pairs))+
 # geom_label_repel(aes(label=Samples, col=pairs),label.size =1, size=5)+
  #geom_point()+
  #labs(x="Proportion of contigs mapped to the genome", y="Total number of contigs (in millions)")+
  #themeAlice+
  #colScale
#myplot2


#Let's remove 2807Do PB CONCENTRATION
newdata <- newdata[-3,]

myplot4<-ggplot(newdata, aes(x=NumMapped, y=PropMapped, label=Samples, col=type, group=type))+
  scale_color_manual(values=c("red","darkgreen","orange"), name="")+
  annotation_custom(grob=circleGrob(r=unit(1,"npc"), gp=gpar(col="red", alpha=1)), xmin=0.02, xmax=0.04, ymin=17.5)+
  geom_point(size=3)+
  geom_label_repel(aes(label=Samples),label.size =1, size=5, show.legend = FALSE)+
  labs(x="Number of contigs mapped to the genome", y="Proportion of contigs mapped to the genome")+
  themeAlice+
  theme(legend.text = element_text(size = 16, face = "bold"),
        legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"),
        legend.position=c(.8, .5))+
  guides(colour = guide_legend(override.aes = list(size=3,linetype=0), title = "Type of procedure"))
myplot4






#*************
#make a distribution plot
#*************

# Histogram overlaid with kernel density curve
Histo <- ggplot(PropDataFrame, aes(x=PropMapped, fill=..count..)) + 
  geom_histogram(binwidth = 0.01, colour="black")+ # Histogram with count on y-axis
  themeAlice+
  scale_x_continuous(name = "Proportion Mapped")+
  scale_y_continuous(breaks = 0:5)+
  scale_fill_continuous(low = "green", high ="red")+
  ggtitle("Distribution of the proportion of mapped sequences among the samples")
Histo
#ggsave("PlotPropMappedDerep.pdf")
#*************


















################################## Later not checked
#1. Total baits assigned
#read.table("/SAN/Alices_sandpit/FeatureCounts_results_dereplicated/OnlyAssigned", sep="\t", header=FALSE)
#length(AssignedBaits$V3)
#   1 273 684 
#Keep only the lines that are duplicated
#HowMany <- duplicated(AssignedBaits$V3)
#table(HowMany)
#Repartition
#  FALSE    TRUE 
#  12391   1261293 
#a <- table(AssignedBaits$V3)
#head(a,10)


#2. Baits assigned per sample
for (i in c(1:8,10,12:14,16:19)) {
  assign(paste("AssignedBaits_",substr(FilesR1[i],1,6), sep = ""), read.table(paste("/SAN/Alices_sandpit/sequencing_data_dereplicated/AssignedFC_",substr(FilesR1[i],1,6), sep = "")))
}

#Find baits similar in several Assignedbaits files

#Compare the third columns of the 19 files
## 1. "WESTERN"
#2808
Inters2808 <- Reduce(intersect, list(AssignedBaits_2808Do$V3,AssignedBaits_2808Si$V3,AssignedBaits_Undete$V3))

#2809
Inters2809 <- Reduce(intersect, list(AssignedBaits_2809Do$V3, AssignedBaits_2809Di$V3)) #AssignedBaits_2809Si$V3 FAILURE
Inters_western_2808_2809 <- intersect(Inters2808, Inters2809) # 567
# put away EfaB_
for (i in 1:length(Inters_western_2808_2809)) {
  Inters_western_2808_2809[i]<-paste(substr(Inters_western_2808_2809[i],6,20))
}

## 2. "EASTERN" Good (<100 contigs mapped)
Inters_eastern_good <- Reduce(intersect, list(AssignedBaits_2672Si$V3, AssignedBaits_2807Di$V3,
                                              AssignedBaits_2807Si$V3, AssignedBaits_2811Do$V3,
                                              AssignedBaits_2812Si$V3, AssignedBaits_2848Si$V3)) # 1
#AssignedBaits_2812Do & AssignedBaits_2919Si FAILURE


# put away EfaB_
for (i in 1:length(Inters_eastern_good)) {
  Inters_eastern_good[i]<-paste(substr(Inters_eastern_good[i],6,20))
}

## 3. "EASTERN" SuperGood (<500 contigs mapped)
Inters_eastern_supergood <- Reduce(intersect, list(AssignedBaits_2672Si$V3, AssignedBaits_2807Di$V3,
                                                   AssignedBaits_2807Si$V3)) #AssignedBaits_2812Do$V3,
                                                   #AssignedBaits_2919Si$V3)) # 31
# put away EfaB_
for (i in 1:length(Inters_eastern_supergood)) {
  Inters_eastern_supergood[i]<-paste(substr(Inters_eastern_supergood[i],6,20))
}

## 4. "EASTERN" TopGood (<500 contigs mapped)
Inters_eastern_topgood <- intersect(AssignedBaits_2672Si$V3, AssignedBaits_2807Si$V3)
Inters_eastern_topgood <- intersect(Inters_eastern_topgood, AssignedBaits_2919Si$V3) # 1948
# put away EfaB_
for (i in 1:length(Inters_eastern_topgood)) {
  Inters_eastern_topgood[i]<-paste(substr(Inters_eastern_topgood[i],6,20))
}

## 5. "apo"
intersApo <- Reduce(intersect, list(AssignedBaits_2TRAnn$V3,AssignedBaits_95Anna$V3,AssignedBaits_f9Anna$V3)) # 116
# put away EfaB_
for (i in 1:length(intersApo)) {
  intersApo[i]<-paste(substr(intersApo[i],6,20))
}

#Find the neigbouring baits!!??


vect <- sort(Inters_western_2808_2809)

for (i in 1:length(vect)) {
  if (substr(vect[1],1,4)==substr(vect[2],1,4))
    listFollowBaits[i]
}
substr(vect[1],1,4)

# Get the bait sequences? [especially the NEIBOURING ONES]




mylist <- lapply(fc2, "[[", "counts") #select fc2$eachsample$counts
mylist2 <- do.call(cbind,mylist) # bind by column
head(mylist2)

a<-table(rowSums(mylist2)>100) #a coverage >100
b<-table(rowSums(mylist2>10)>5) #a coverage >10 in >5 lib
c<-table(rowSums(mylist2>0)>15) #a coverage >0 in >15 lib
names(dimnames(a)) <- list("coverage >100")
names(dimnames(b)) <- list("coverage >10 in >5 lib")
names(dimnames(c)) <- list("coverage >0 in >15 lib")
tot <- cbind(a,b,c)
colnames(tot)<- c("coverage >100","coverage >10 in >5 lib","coverage >0 in >15 lib")
tot
tot<-xtable(tot, align = c("c","c","c","c"))
print.xtable(tot, file="tabcoverage.html", type = "html")





#**Citation**
#  Yang Liao, Gordon K Smyth and Wei Shi (2013). The Subread aligner: fast, accurate
#and scalable read mapping by seed-and-vote. Nucleic Acids Research, 41(10):e108.
#Yang Liao, Gordon K Smyth and Wei Shi (2014). featureCounts: an efficient general
#purpose program for assigning sequence reads to genomic features. Bioinformatics,
#30(7):923-30

#** RE general **
#  http://bioinformatics.oxfordjournals.org/content/30/7/923.full.pdf?keytype=ref&ijkey=ZzPz96t2lqzAH6F
#http://bioconductor.org/packages/3.3/bioc/vignettes/ShortRead/inst/doc/Overview.pdf
#https://darrenjw.wordpress.com/tag/fastq/
#  To install further bioconductor packages :
#  source("http://bioconductor.org/biocLite.R")  
#biocLite("flowViz") 
#http://www.mycroarray.com/mybaits/mybaits-technology.html