###libraries
## setwd("/SAN/Alices_sandpit/sequencing_data_dereplicated/")
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
library(reshape)
library(viridis)
library(ggthemes)
library(pheatmap)
library(GGally)

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


#####  I

            
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
#Tallying implicitly paired files
#In this scenario two files are implicitly paired, such as is the case for unprocessed paired-end FASTQ files. Tally will process both files record-by-record and pair up records at the same offset. This requires the option --pair-by-offset.
#tally -i out1.gz -j out2.gz -o out1.unique.gz -p out2.unique.gz --pair-by-offset

#And parallel...... to be continued

#Rsubread package: (3) align the dereplicated sequences to the index
#Map paired-end reads
FilesR1 <- list.files(path = "/SAN/Alices_sandpit/sequencing_data_dereplicated",
                      pattern="R1_001.fastq.unique.gz", full.names=TRUE)
FilesR2 <- list.files(path = "/SAN/Alices_sandpit/sequencing_data_dereplicated",
                      pattern="R2_001.fastq.unique.gz", full.names=TRUE)
FilesR1; FilesR2

#Played around with snow, really tough, ask Totta for help with parallelisation proccess

## Run again with slimmer Code ## try using mclapply (multi-core apply) or snow for parallelization...
## library(snow) or library(mcapply)

#Create 1 alignments (20min/alignment !!)
#*************
####for (i in 1:length(FilesR1path)){
####  reads1 <- FilesR1path[i]
####  reads2 <- FilesR2path[i]
####  Rsubread::align(index="reference_index",readfile1=reads1,readfile2=reads2, type="dna",maxMismatches=20, indels=10,
####                  output_file = paste("Alignment_", substr(FilesR1[i],1,6), sep = ""))
####}

# Add a . bam !!!!!

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



#Cf Totta things: 
#parallel --gnu -P 5 --xapply tophat -r 200 --library-type fr-unstranded -G reference_genomes/indexes_bowtie2_good/mm10_GRCm38_eimeriaHaberkorn.gtf -o tophat_March_c/{1/.}_paired reference_genomes/indexes_bowtie2_good/index_mm10_eimeria3 {1} {2} ::: /data/Eimeria_Totta/RNAseq*/*_forw.fastq.gz ::: /data/Eimeria_Totta/RNAseq*/*_rev.fastq.gz


######################
##Mapping percentage##
######################
#Function propmapped returns the proportion of mapped reads included in a SAM/BAM
#file. For paired end reads, it can return the proportion of mapped fragments (ie. read pairs).
#*************

Align.file <- list.files(path = "/SAN/Alices_sandpit/sequencing_data_dereplicated",
                         pattern="Alignment_", full.names=TRUE)

Prop.mapped <- propmapped(Align.file)

Prop.mapped$Samples <- substr(Prop.mapped$Samples,60,65)
Prop.mapped$propMapped <- Prop.mapped$PropMapped*100

# Histogram overlaid with kernel density curve
Histo <- ggplot(Prop.mapped, aes(x=propMapped)) + 
  geom_histogram(colour="black", fill="white")+ # Histogram with count on y-axis
  themeAlice+
  scale_x_continuous(name = "Proportion Mapped")+
  scale_y_continuous(breaks = 0:10)
Histo

###############################################
## Counting mapped reads for genomic features##
###############################################
fc2 <- list()
FilesR1
for (i in  c(1:8,10,12:14,16:19)) {
  fc2[[FilesR1[i]]] <-  featureCounts(paste("/SAN/Alices_sandpit/sequencing_data_dereplicated/Alignment_", substr(FilesR1[i],50,55), sep = ""),
                                      annot.ext="/SAN/Alices_sandpit/sequencing_data_dereplicated/MYbaits_Eimeria_V1.single120_feature_counted.gtf",
                                      isGTFAnnotationFile=TRUE,
                                      GTF.featureType = "sequence_feature",
                                      useMetaFeatures=FALSE,
                                      GTF.attrType="bait_id",
                                      # parameters specific to paired end reads
                                      isPairedEnd=TRUE,
                                      reportReads=TRUE)
}

# First look into the data
mylist <- lapply(fc2, "[[", "counts")
countsDF <- do.call(cbind,mylist)
cgountsDF <- as.data.frame(countsDF)
colnames(countsDF) <-  substr(colnames(countsDF),61,66) # shorten names 
table(rowSums(countsDF)>100)
table(rowSums(countsDF>10)>5) #a coverage >10 in >5 lib = how many baits are "working" with more than 10 sequences captured in more than 5 libraries
table(rowSums(countsDF>0)>10)#a coverage >0 in >10 lib = how many baits are "working" in more than 10 libraries

#For later: change the names of countsDF libraries with a "Lib" in front
colnames(countsDF)<- paste("Lib_",colnames(countsDF),sep="")

###################################################
#Hierarchical clustering on the readcounts per bait.
#A heatmap just on this raw data?
#Pairs plot to identify correlating (good) libraries.
#Table of pairwise correlation coefficients.
#And a truly random 120nt baits gff would be awesome ;-)
#Just random no intron or etc selection needed.
#My guess is the hierarchical clustering will show a cluster of "working" baits. And the correlations will show "working" libraries.
#column scaling of the libraries might give the best clustering.
################################################### 

#Plot a heatmap
pheatmap(log10(countsDF[rowSums(countsDF)>500,]+0.1))
pheatmap(log10(countsDF[rowSums(countsDF)>50,]+0.1))
pheatmap(log10(countsDF[rowSums(countsDF)>10,]+0.1))
pheatmap(log10(countsDF[rowSums(countsDF)>500,]+0.1))
table(rowSums(countsDF)>10)

############
###GOOD library, definition
############
minbaits <- 5000 # we want minimum b "good baits"
mincounts <- 10 # with minimum c counts each
totbaits <- nrow(countsDF)
############

# Keep the lib with at least b baits with greater than c counts
goodlib <- apply(countsDF,2, function(x) sum((table(x))[c(1:mincounts)])<totbaits-minbaits) #if NA, there is always less than c counts per bait/if FALSE there is not enough "good" baits/we want to keep the TRUE ones
NewDF <- countsDF[which(goodlib==TRUE)]
#Now, remove rows with zeros on the lib selected
##Go through each row and determine if a value is zero
row_sub <- apply(NewDF, 1, function(row) all(row !=0 ))
##Subset as usual
FinalDF <- NewDF[row_sub,]
baits.to.keep <- nrow(FinalDF)
head(FinalDF)
#################
##### countsDF is the original table,
##### NewDF is the table with just "good" libraries
##### FinalDF is the table with "good" libraries and no zeros!
#################

#Pairs plot to identify correlating (good) libraries
my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) +
    geom_point() +
    geom_smooth(method=loess, fill="red", color="red", ...) +
    geom_smooth(method=lm, fill="blue", color="blue", ...)
  p
}

#pdf("/SAN/Alices_sandpit/sequencing_data_dereplicated/figures/countsDF_ggpairs.pdf", height = 7, width = 7)
#g <- ggpairs(countsDF, lower = list(continuous = my_fn))
#print(g)
#dev.off()
                                        # Saved in Alice_sandpit/sequencing_data_dereplicated/figures


###############################
# Plot nice correlation heatmap + pairwise coef
#Function adapted from http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
###############################
create_heatmap_correl <- function(mydata){
                                        # I.Compute the correlation matrix
cormat <- round(cor(mydata),2)
                                        # II.Get upper triangle of the correlation matrix
    get_upper_tri <- function(cormat){
      cormat[lower.tri(cormat)]<- NA
      return(cormat)
    }
  upper_tri <- get_upper_tri(cormat)
                                        # II.Melt the correlation matrix
  melted_cormat <- na.omit(melt(upper_tri))
                                        # IV.Plot
    ggplot(data = melted_cormat, aes(X2,X1, fill = value))+
      geom_tile(color = "white")+
      scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                                 midpoint = 0, limit = c(-1,1), space = "Lab",
                                 name="Pearson\nCorrelation") +
      theme_minimal()+
      theme(axis.text.x = element_text(angle = 45, vjust = 1,
              size = 12, hjust = 1))+
                coord_fixed()+                                        #Add correlation coefficients on the heatmap
      geom_text(aes(X2, X1, label = value), color = "black", size = 4) +
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.6, 0.7),
        legend.direction = "horizontal")+
      guides(fill = guide_colorbar(barwidth = 7, barheight = 1,title.position = "top", title.hjust = 0.5))
             }

create_heatmap_correl(countsDF)
create_heatmap_correl(NewDF)
create_heatmap_correl(FinalDF)

                                        # On create_heatmap_correl(countsDF) we can identify "bad" libraries
##Enter here the "bad" libraries / define a criterium for it: here the 3 worst
countsDFreduced <- countsDF[ , -which(names(countsDF) %in% c("Lib_2807Do","Lib_2919Di","Lib_f9Anna"))]

create_heatmap_correl(countsDFreduced)

###################################################
#Hierarchical clustering on the readcounts per bait.  V
#A heatmap just on this raw data?  V
#Pairs plot to identify correlating (good) libraries.  V
#Table of pairwise correlation coefficients.  V
#And a truly random 120nt baits gff would be awesome ;-)
#Just random no intron or etc selection needed.
#My guess is the hierarchical clustering will show a cluster of "working" baits. And the correlations will show "working" libraries.
#column scaling of the libraries might give the best clustering.
################################################### 
                                        #is there reads to non baits region in the genomes?















###### Plots & Stuff Not optimized DO NOT RUN 
#
###### PART II
#PropDataFrame <- propmapped(paste("Alignment_", substr(FilesR1[1],1,6), sep = ""))
#for (i in c(2:8,10,12:14,16:19)) {
#  PropDataFrame <- rbind(PropDataFrame, propmapped(paste("Alignment_", substr(FilesR1[i],1,6), sep = "")))
#}
#
## Number of mapped reads/fragments will be counted and fraction of such reads/fragments will be calculated.
#
#
## Save as PDF
#PropDataFrame
#library(gridExtra)
#pdf("PropDataFrame.pdf", height=11, width=8.5)
#grid.table(PropDataFrame)
#dev.off()
##
#
## :) Shorten the names for better visualisation
#PropDataFrame$Samples<-base::gsub("Alignment_","",PropDataFrame$Samples)
#
#newdata <- PropDataFrame[base::order(PropDataFrame$Samples),]
##add a column to see if single, double, digested
#for (i in 1:nrow(newdata)){
#  if (grepl("Do",newdata$Samples[i])) {
#    newdata$type[i] <- "Double"
#  }
#  else if (grepl("Si",newdata$Samples[i])) {
#    newdata$type[i] <- "Single"
#  }
#  else if (grepl("Di",newdata$Samples[i])) {
#    newdata$type[i] <- "Digested"
#  }
#}
#
##!!! "undetermine" is "2808Di"
#newdata$Samples[16] <- "2808Di"
#newdata$type[16] <- "Digested"
#
## Nbr in million
#newdata$NumTotalinmillion <-newdata$NumTotal/(1000000)
#
#myplot<-ggplot(newdata, aes(x=PropMapped, y=NumTotalinmillion, label=Samples, col=type))+
#  scale_color_manual(values=c(3,2,1), name="")+
#  geom_point()+
#  geom_label_repel(aes(label=Samples),label.size =1, size=5)+
#  labs(x="Proportion of contigs mapped to the genome", y="Total number of contigs (in millions)")+
#  themeAlice
##annotation_custom(grob=circleGrob(r=unit(1,"npc"), gp=gpar(fill="orange", alpha=0.2)), xmin=0.161, xmax=0.179, ymin=10.1, ymax=15.9)+
##annotate("text", x = rep(0.17,3), y = 12:14, size=8, label = c("Digested","Double","Single"),col=c(3,2,1))
#myplot
#
############################################
##add a column to see Eastern/Western clade#
############################################
##2672 Eastern musculus
##2807 Eastern domesticus
##2808 Western domesticus
##2809 Western domesticus
##2811 Eastern domesticus
##2812 Eastern domesticus
##2848 Eastern hybrid
##2919 bEastern domesticus
##2TRAnn X apodemus
##95Anna X apodemus
##f9Anna X apodemus
#
## East/west classif: 
#newdata$clade <- ifelse(grepl("2672",newdata$Samples),"Eastern",NA)
#newdata$clade <- ifelse(grepl("2807",newdata$Samples),"Eastern",newdata$clade)
#newdata$clade <- ifelse(grepl("2808",newdata$Samples),"Western",newdata$clade)
#newdata$clade <- ifelse(grepl("2809",newdata$Samples),"Western",newdata$clade)
#newdata$clade <- ifelse(grepl("2811",newdata$Samples),"Eastern",newdata$clade)
#newdata$clade <- ifelse(grepl("2812",newdata$Samples),"Eastern",newdata$clade)
#newdata$clade <- ifelse(grepl("2848",newdata$Samples),"Eastern",newdata$clade)
#newdata$clade <- ifelse(grepl("2919",newdata$Samples),"Eastern",newdata$clade)
#newdata$clade <- ifelse(grepl("2TRAnn",newdata$Samples),"X",newdata$clade)
#newdata$clade <- ifelse(grepl("95Anna",newdata$Samples),"X",newdata$clade)
#newdata$clade <- ifelse(grepl("f9Anna",newdata$Samples),"X",newdata$clade)
#
## Mus/dom/Apo classif: 
#newdata$mice <- ifelse(grepl("2672",newdata$Samples),"musculus",NA)
#newdata$mice <- ifelse(grepl("2807",newdata$Samples),"domesticus",newdata$mice)
#newdata$mice <- ifelse(grepl("2808",newdata$Samples),"domesticus",newdata$mice)
#newdata$mice <- ifelse(grepl("2809",newdata$Samples),"domesticus",newdata$mice)
#newdata$mice <- ifelse(grepl("2811",newdata$Samples),"domesticus",newdata$mice)
#newdata$mice <- ifelse(grepl("2812",newdata$Samples),"domesticus",newdata$mice)
#newdata$mice <- ifelse(grepl("2848",newdata$Samples),"hybrid",newdata$mice)
#newdata$mice <- ifelse(grepl("2919",newdata$Samples),"domesticus",newdata$mice)
#newdata$mice <- ifelse(grepl("2TRAnn",newdata$Samples),"apodemus",newdata$mice)
#newdata$mice <- ifelse(grepl("95Anna",newdata$Samples),"apodemus",newdata$mice)
#newdata$mice <- ifelse(grepl("f9Anna",newdata$Samples),"apodemus",newdata$mice)
#
#
#myplot3<-ggplot(newdata, aes(x=PropMapped, y=NumTotalinmillion, label=Samples, col=clade))+
#  scale_color_manual(values=c("green","darkgreen","orange"), name="")+
#  annotation_custom(grob=circleGrob(r=unit(1,"npc"), gp=gpar(col="red", alpha=1)), xmin=0.02, xmax=0.04, ymin=17.5)+
#  geom_point(size=3)+
#  geom_label_repel(aes(label=Samples),label.size =1, size=5)+
#  labs(x="Proportion of contigs mapped to the genome", y="Total number of contigs (in millions)")+
#  themeAlice+
#  annotate("text", x = rep(0.16,3), y = 12:14, size=8, label = c("Closely related to Western clade","Eastern clade","Western clade"),col=c("green","darkgreen","orange"))
#myplot3
#
####### ggplot function
##Create a custom color scale
#myColors <- brewer.pal(12,"Paired")
#names(myColors) <- levels(newdata$pairs)
#colScale <- scale_colour_manual(name = "pairs",values = myColors)
#######
#
##add a column to see the pairs
##newdata$pairs<- c(1,2,2,2,3,3,4,4,4,5,6,6,7,8,8,9,10,11,12)
##newdata$pairs <- as.factor(newdata$pairs)
#
##Plot
##myplot2<-ggplot(newdata, aes(x=PropMapped, y=NumTotalinmillion, label=Samples, col=pairs))+
# # geom_label_repel(aes(label=Samples, col=pairs),label.size =1, size=5)+
#  #geom_point()+
#  #labs(x="Proportion of contigs mapped to the genome", y="Total number of contigs (in millions)")+
#  #themeAlice+
#  #colScale
##myplot2
#
##Let's remove 2807Do PB CONCENTRATION
#newdata <- newdata[-3,]
#
#myplot4<-ggplot(newdata, aes(x=NumMapped, y=PropMapped, label=Samples, col=type, group=type))+
#  scale_color_manual(values=c("red","darkgreen","orange"), name="")+
#  annotation_custom(grob=circleGrob(r=unit(1,"npc"), gp=gpar(col="red", alpha=1)), xmin=0.02, xmax=0.04, ymin=17.5)+
#  geom_point(size=3)+
#  geom_label_repel(aes(label=Samples),label.size =1, size=5, show.legend = FALSE)+
#  labs(x="Number of contigs mapped to the genome", y="Proportion of contigs mapped to the genome")+
#  themeAlice+
#  theme(legend.text = element_text(size = 16, face = "bold"),
#        legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"),
#        legend.position=c(.8, .5))+
#  guides(colour = guide_legend(override.aes = list(size=3,linetype=0), title = "Type of procedure"))
#myplot4
#
##*************
##make a distribution plot
##*************
#
## Histogram overlaid with kernel density curve
#Histo <- ggplot(PropDataFrame, aes(x=PropMapped, fill=..count..)) + 
#  geom_histogram(binwidth = 0.01, colour="black")+ # Histogram with count on y-axis
#  themeAlice+
#  scale_x_continuous(name = "Proportion Mapped")+
#  scale_y_continuous(breaks = 0:5)+
#  scale_fill_continuous(low = "green", high ="red")+
#  ggtitle("Distribution of the proportion of mapped sequences among the samples")
#Histo
##ggsave("PlotPropMappedDerep.pdf")
##*************
#
#1. Total baits assigned
##read.table("/SAN/Alices_sandpit/FeatureCounts_results_dereplicated/OnlyAssigned", sep="\t", header=FALSE)
##length(AssignedBaits$V3)
##   1 273 684 
##Keep only the lines that are duplicated
##HowMany <- duplicated(AssignedBaits$V3)
##table(HowMany)
##Repartition
##  FALSE    TRUE 
##  12391   1261293 
##a <- table(AssignedBaits$V3)
##head(a,10)
#
#
##2. Baits assigned per sample
#for (i in c(1:8,10,12:14,16:19)) {
#  assign(paste("AssignedBaits_",substr(FilesR1[i],1,6), sep = ""), read.table(paste("/SAN/Alices_sandpit/sequencing_data_dereplicated/AssignedFC_",substr(FilesR1[i],1,6), sep = "")))
#}
#
##Find baits similar in several Assignedbaits files
#
##Compare the third columns of the 19 files
### 1. "WESTERN"
##2808
#Inters2808 <- Reduce(intersect, list(AssignedBaits_2808Do$V3,AssignedBaits_2808Si$V3,AssignedBaits_Undete$V3))
#
##2809
#Inters2809 <- Reduce(intersect, list(AssignedBaits_2809Do$V3, AssignedBaits_2809Di$V3)) #AssignedBaits_2809Si$V3 FAILURE
#Inters_western_2808_2809 <- intersect(Inters2808, Inters2809) # 567
## put away EfaB_
#for (i in 1:length(Inters_western_2808_2809)) {
#  Inters_western_2808_2809[i]<-paste(substr(Inters_western_2808_2809[i],6,20))
#}
#
### 2. "EASTERN" Good (<100 contigs mapped)
#Inters_eastern_good <- Reduce(intersect, list(AssignedBaits_2672Si$V3, AssignedBaits_2807Di$V3,
#                                              AssignedBaits_2807Si$V3, AssignedBaits_2811Do$V3,
#                                              AssignedBaits_2812Si$V3, AssignedBaits_2848Si$V3)) # 1
##AssignedBaits_2812Do & AssignedBaits_2919Si FAILURE
#
#
## put away EfaB_
#for (i in 1:length(Inters_eastern_good)) {
#  Inters_eastern_good[i]<-paste(substr(Inters_eastern_good[i],6,20))
#}
#
### 3. "EASTERN" SuperGood (<500 contigs mapped)
#Inters_eastern_supergood <- Reduce(intersect, list(AssignedBaits_2672Si$V3, AssignedBaits_2807Di$V3,
#                                                   AssignedBaits_2807Si$V3)) #AssignedBaits_2812Do$V3,
#                                                   #AssignedBaits_2919Si$V3)) # 31
## put away EfaB_
#for (i in 1:length(Inters_eastern_supergood)) {
#  Inters_eastern_supergood[i]<-paste(substr(Inters_eastern_supergood[i],6,20))
#}
#
### 4. "EASTERN" TopGood (<500 contigs mapped)
#Inters_eastern_topgood <- intersect(AssignedBaits_2672Si$V3, AssignedBaits_2807Si$V3)
#Inters_eastern_topgood <- intersect(Inters_eastern_topgood, AssignedBaits_2919Si$V3) # 1948
## put away EfaB_
#for (i in 1:length(Inters_eastern_topgood)) {
#  Inters_eastern_topgood[i]<-paste(substr(Inters_eastern_topgood[i],6,20))
#}
#
### 5. "apo"
#intersApo <- Reduce(intersect, list(AssignedBaits_2TRAnn$V3,AssignedBaits_95Anna$V3,AssignedBaits_f9Anna$V3)) # 116
## put away EfaB_
#for (i in 1:length(intersApo)) {
#  intersApo[i]<-paste(substr(intersApo[i],6,20))
#}
#
##Find the neigbouring baits!!??
#
#
#vect <- sort(Inters_western_2808_2809)
#
#for (i in 1:length(vect)) {
#  if (substr(vect[1],1,4)==substr(vect[2],1,4))
#    listFollowBaits[i]
#}
#substr(vect[1],1,4)
#
## Get the bait sequences? [especially the NEIBOURING ONES]
#
#
#a<-table(rowSums(mylist2)>100) #a coverage >100
#b<-table(rowSums(mylist2>10)>5) #a coverage >10 in >5 lib
#c<-table(rowSums(mylist2>0)>15) #a coverage >0 in >15 lib
#names(dimnames(a)) <- list("coverage >100")
#names(dimnames(b)) <- list("coverage >10 in >5 lib")
#names(dimnames(c)) <- list("coverage >0 in >15 lib")
#tot <- cbind(a,b,c)
#colnames(tot)<- c("coverage >100","coverage >10 in >5 lib","coverage >0 in >15 lib")
#tot
#tot<-xtable(tot, align = c("c","c","c","c"))
#print.xtable(tot, file="tabcoverage.html", type = "html")


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
