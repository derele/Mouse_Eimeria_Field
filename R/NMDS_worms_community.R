## NMDS on mice and their worms community
## Alice 21 November 2017
## Tuto in https://jonlefcheck.net/2012/10/24/nmds-tutorial-in-r/

source("Data_Analysis_Alice_2017.R")

# Non-metric multidimensional scaling (NMDS) is one tool commonly used to 
# examine community composition

# install.packages("vegan")
library(vegan)
# `metaMDS` requires a community-by-species matrix
worms <- TotalTable[c("Cysticercus", "Trichuris_muris", "Aspiculuris_tetraptera", "Syphacia_obvelata",
                      "Mastophorus_muris", "Heterakis_spumosa", "Mesocestoides", 
                      "Catenotaenia_pusilla", "Hymenolepis", "Oxyurids", "Mix_Syphacia_Aspiculuris",
                      "Heligmosomoides_polygurus", "HI", "Mouse_ID", "mice.per.trap")]
# ("Hymenolepis_microstoma", "Hymenolepis_diminiuta") to check, contains "TRUE" and "FALSE"

# remove non mus musculus
worms <- worms[-grep("ZZ_", worms$Mouse_ID),]

### 1. organise by subspecies
worms$subspe <- NA
worms$subspe[!is.na(worms$HI)] <- "hybrids"
worms$subspe[!is.na(worms$HI) & (worms$HI < 0.1)]  <- "Mmd"
worms$subspe[!is.na(worms$HI) & (worms$HI > 0.9)]  <- "Mmm"

my_community_matrix <- data.matrix(worms)
dimnames(my_community_matrix)[1] <- list(worms[,"subspe"])

# remove rows with only NA for worms (not checked for worms)
my_community_matrix <-
  my_community_matrix[rowSums(is.na(my_community_matrix[,c(1:12)]))!= ncol(my_community_matrix[,c(1:12)]), ]

# keep worms only
my_community_matrix <- my_community_matrix[,c(1:12)]

# then add zeros when NA for worms (absence of worm species checked)
my_community_matrix[is.na(my_community_matrix)] <- 0

# data normalization : (x - mean)/sd
my_community_matrix <- apply(my_community_matrix, 2, function(x) {(x - mean(x))/ sd(x)})

# avoid the negative data problem
my_community_matrix <- my_community_matrix + 1

# The function `metaMDS` will take care of most of the distance 
# calculations, iterative fitting, etc. We need simply to supply:
my_NMDS = metaMDS(comm = my_community_matrix, k = 2, maxit = 1000, sratmax = 100)
##**A good rule of thumb: stress < 0.05 provides an excellent representation in reduced dimensions,
## < 0.1 is great, < 0.2 is good/ok, and stress < 0.3 provides a poor representation.** 

# Shepard plot, which shows scatter around the regression
# between the interpoint distances in the final configuration (distances 
# between each pair of communities) against their original dissimilarities
stressplot(my_NMDS)
# Large scatter around the line suggests that original dissimilarities are
# not well preserved in the reduced number of dimensions

#Now we can plot the NMDS
plot(my_NMDS)

col = rep("black",dim(my_community_matrix)[1])
r = rownames(my_community_matrix)
col[r=="hybrids"] = "violet"
col[r=="Mmm"] = "red"
col[r=="Mmd"] = "blue"

ordiplot(my_NMDS,type="n")
orditorp(my_NMDS,display="species",col="black",air=0.01)
orditorp(my_NMDS,display="sites",col=col,air=0.01,cex=1.25)

#######################################
## 2. organise by mice.per.traps groups
hist(TotalTable$mice.per.trap, breaks = 50)
# 2 groups : high density, low density (<40, >40)

# `metaMDS` requires a community-by-species matrix
worms <- TotalTable[c("Cysticercus", "Trichuris_muris", "Aspiculuris_tetraptera", "Syphacia_obvelata",
                      "Mastophorus_muris", "Heterakis_spumosa", "Mesocestoides", 
                      "Catenotaenia_pusilla", "Hymenolepis", "Oxyurids", "Mix_Syphacia_Aspiculuris",
                      "Heligmosomoides_polygurus", "HI", "Mouse_ID", "mice.per.trap")]
# ("Hymenolepis_microstoma", "Hymenolepis_diminiuta") to check, contains "TRUE" and "FALSE"

# remove non mus musculus
worms <- worms[-grep("ZZ_", worms$Mouse_ID),]

worms$mice.dens <- " "
worms$mice.dens[worms$mice.per.trap < 40]  <- "low"
worms$mice.dens[worms$mice.per.trap > 40]  <- "high"

my_community_matrix <- data.matrix(worms)
dimnames(my_community_matrix)[1] <- list(worms[,"mice.dens"])

# remove rows with only NA for worms (not checked for worms)
my_community_matrix <-
  my_community_matrix[rowSums(is.na(my_community_matrix[,c(1:12)]))!= ncol(my_community_matrix[,c(1:12)]), ]

# keep worms only
my_community_matrix <- my_community_matrix[,c(1:12)]

# then add zeros when NA for worms (absence of worm species checked)
my_community_matrix[is.na(my_community_matrix)] <- 0

# data normalization : (x - mean)/sd
my_community_matrix <- apply(my_community_matrix, 2, function(x) {(x - mean(x))/ sd(x)})

# avoid the negative data problem
my_community_matrix <- my_community_matrix + 1

# The function `metaMDS` will take care of most of the distance 
# calculations, iterative fitting, etc. We need simply to supply:
my_NMDS = metaMDS(comm = my_community_matrix, k = 2, maxit = 1000, sratmax = 100)
##**A good rule of thumb: stress < 0.05 provides an excellent representation in reduced dimensions,
## < 0.1 is great, < 0.2 is good/ok, and stress < 0.3 provides a poor representation.** 

# Shepard plot, which shows scatter around the regression
# between the interpoint distances in the final configuration (distances 
# between each pair of communities) against their original dissimilarities
stressplot(my_NMDS)
# Large scatter around the line suggests that original dissimilarities are
# not well preserved in the reduced number of dimensions

#Now we can plot the NMDS
plot(my_NMDS)

col = rep("black",dim(my_community_matrix)[1])
r = rownames(my_community_matrix)
col[r==" "] = "lightgrey"
col[r=="low"] = "skyblue"
col[r=="high"] = "darkred"

ordiplot(my_NMDS,type="n")
orditorp(my_NMDS,display="species",col="black",air=0.01)
orditorp(my_NMDS,display="sites",col=col,air=0.01,cex=1.25)