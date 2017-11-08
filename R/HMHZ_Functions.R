## list of usefull function August 2017 :
## E. Heitlinger
## A. Balard

library(ggmap)

## Convert GPS coordinate in decimals format
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

## Calculate the hybrid index
get.HIX <- function (x){
    dom <- nchar(x) - nchar(gsub("d", "", x))
    mus <- nchar(x) - nchar(gsub("m", "", x))
    mus/(mus + dom)
}

is.loc.cluster <- function(locA, locB, granularity=0.001){
    diff <- as.numeric(as.character(locA[c("Latitude", "Longitude")])) -
        as.numeric(as.character(locB[c("Latitude", "Longitude")]))
    ## both east-west and north-south are closer than granularit
    all (abs(diff) < granularity)
}

## pairwise comparison of all localities 
pairwise.cluster.loc <- function (d){
    loc.combis <- combn(nrow(d), 2)
    long <- apply(loc.combis, 2, function(x)
        is.loc.cluster(d[x[1], ], d[x[2], ]))
    mat <- matrix(NA, nrow=nrow(d), ncol=nrow(d))
    mat[lower.tri(mat)] <- long
    mat <- t(mat)
    ## setting lower triangle and diagonal zero
    diag(mat) <- 0
    mat[lower.tri(mat)] <- 0
    mat
}

foo <- pairwise.cluster.loc(Trap)

## how many locs are sampled multiple times
table(apply(foo, 2, sum))

## Everyting below here is gibberish and should be rewritten ;-)

## map rows from different data.frames by their common column
rbind.match.columns <- function(input1, input2) {
    n.input1 <- ncol(input1)
    n.input2 <- ncol(input2)
    if (n.input2 < n.input1) {
        TF.names <- which(names(input2) %in% names(input1))
        column.names <- names(input2[, TF.names])
    } else {
        TF.names <- which(names(input1) %in% names(input2))
        column.names <- names(input1[, TF.names])
    }
    return(rbind(input1[, column.names], input2[, column.names]))
}

## Create a genotype dataframe :
make.gen.DF <- function(){
  
  # Add the last data received :
  genotypes.2016.and.some.previous <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/raw_data/HIforEH_May2017.csv") 

  Gen.almost.tot <- data.frame(Mouse_ID = genotypes.2016.and.some.previous$PIN, Code = genotypes.2016.and.some.previous$Code,
                               Year = genotypes.2016.and.some.previous$Year, HI = genotypes.2016.and.some.previous$HI,
                               HI_NLoci = genotypes.2016.and.some.previous$HI_NLoci)
  Gen.almost.tot$Mouse_ID <- gsub(pattern = "SK", replacement = "SK_",x = Gen.almost.tot$Mouse_ID)
  
  # Previous data
  genotypes.2014 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/raw_data/HZ14_Mice%2031-12-14_genotypes.csv")[-1,]
  X.col <- c("X332", "X347", "X65", "Tsx", "Btk",  "Syap1")
  genotypes.2014$Mouse_ID <- paste(genotypes.2014$ID,genotypes.2014$PIN, sep = "_")
  
  genotypes.2015.1 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/raw_data/HZ15_Mice_genotypes.csv")
  genotypes.2015.1$Mouse_ID <- paste(genotypes.2015.1$ID,genotypes.2015.1$PIN, sep = "_")
  genotypes.2015.1$Year <- "2015"
  
  genotypes.2015.2 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/raw_data/Genotypes_Bav2015.csv")
  genotypes.2015.2$Mouse_ID <- gsub(pattern = "SK", replacement = "SK_",x = genotypes.2015.2$PIN)
  
  # Bind them all :
  MyGen <- function(x) {
    cbind(x[,X.col], x$Mouse_ID, x$Code, x$Year)
  }
  
  Gen14.15 <- rbind(MyGen(genotypes.2014),
                    MyGen(genotypes.2015.1),
                    MyGen(genotypes.2015.2))
  
  ## collapse m or d alleles in 1 column
  Gen14.15$collapsed.GT <-
    apply(Gen14.15[, X.col],1, paste, collapse = "/")
  
  ## calculate hybrid index per mouse
  Gen14.15$HIX <-
    sapply(Gen14.15$collapsed.GT, get.HIX)
  
  Gen14.15 <- Gen14.15[c(7, 8, 9, 11)] 
  Gen14.15$HI_NLoci <- "HI 6"
  names(Gen14.15) <- c("Mouse_ID", "Code", "Year", "HI", "HI_NLoci")
  
  # Manual correction 
  Gen14.15[grep("Sk", Gen14.15$Mouse_ID), ]$Mouse_ID <- "SK_3173"
 
  # Merge both data frames
  Gen.tot <- merge(Gen.almost.tot, Gen14.15, by = c("Mouse_ID", "Code", "Year"), all = TRUE)
  
  # By default :
  Gen.tot$HI <- NA
  Gen.tot$HI_NLoci <- NA
  Gen.tot[which(!is.na(Gen.tot$HI.y)),]$HI <- Gen.tot[which(!is.na(Gen.tot$HI.y)),]$HI.y
  Gen.tot[which(!is.na(Gen.tot$HI_NLoci.y)),]$HI_NLoci <- Gen.tot[which(!is.na(Gen.tot$HI_NLoci.y)),]$HI_NLoci.y
  
  # But if "better" HI calculated :
  Gen.tot[which(!is.na(Gen.tot$HI.x)),]$HI <- Gen.tot[which(!is.na(Gen.tot$HI.x)),]$HI.x
  Gen.tot[which(!is.na(Gen.tot$HI_NLoci.x)),]$HI_NLoci <- as.character(Gen.tot[which(!is.na(Gen.tot$HI_NLoci.x)),]$HI_NLoci.x)
  
  # Keep best HI if one calculated with more markers
  Gen.tot$HI_NLoci <- as.numeric(gsub(pattern = "HI ", replacement = "", x = Gen.tot$HI_NLoci))
  
  # Remove useless 
  Gen.tot <- Gen.tot[-(4:7)]
  
  Gen.tot <- Gen.tot[order(Gen.tot$Mouse_ID, -Gen.tot$HI), ] #sort by id and reverse of HI
  Gen.tot <- Gen.tot[!duplicated(Gen.tot$Mouse_ID), ]              # take the first row within each id
  
  Gen.tot <- na.omit(Gen.tot)
  Gen.tot <- unique(Gen.tot)
  
  # Last check duplicates
  Gen.tot[which(duplicated(Gen.tot$Mouse_ID)),]

  # Add latitude / longitude
  
  # All latitude / longitude unknown should be in the file Gen.almost.tot
  Jaro_loc <- unique(rbind(unique(data.frame(Code = genotypes.2016.and.some.previous$Code, 
                                           Longitude = genotypes.2016.and.some.previous$Xmap, 
                                          Latitude = genotypes.2016.and.some.previous$Ymap))))
  
  # Cleaned locations
  cleaned_loc <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/output_data/all_clean_localities.csv")
  cleaned_loc <- cleaned_loc[-1]
  cleaned_loc$Code <- gsub(pattern = "E_", replacement = "", x = cleaned_loc$Code)
  
  # Merge
  All_loc <- merge(cleaned_loc, Jaro_loc, by = "Code", all = TRUE)
  
  # By default take cleaned values 
  All_loc$Longitude <- NA
  All_loc$Latitude <- NA
  All_loc[which(!is.na(All_loc$Latitude.x)),]$Latitude <- All_loc[which(!is.na(All_loc$Latitude.x)),]$Latitude.x 
  All_loc[which(!is.na(All_loc$Longitude.x)),]$Longitude <- All_loc[which(!is.na(All_loc$Longitude.x)),]$Longitude.x 
  
  # Otherwise take Jaroslav coordinates
  All_loc[which(is.na(All_loc$Latitude)),]$Latitude <- All_loc[which(is.na(All_loc$Latitude)),]$Latitude.y
  All_loc[which(is.na(All_loc$Longitude)),]$Longitude <- All_loc[which(is.na(All_loc$Longitude)),]$Longitude.y 
  
  # Give rank
  All_loc$rank <- 2
  All_loc[which(is.na(All_loc$Latitude.x)),]$rank <- 2
  
  # Delete duplicates based on rank
  All_loc <- All_loc[order(All_loc$Code, -All_loc$rank), ] #sort by id and reverse of HI
  All_loc <- All_loc[!duplicated(All_loc$Code), ]              # take the first row within each id
  All_loc <- All_loc[-c(2:5,8)]
  
  # Finally, merge Gen and Loc
  Gen.and.loc <- merge(Gen.tot, All_loc, by = c("Code"), all = TRUE)
  
  # For which code and year do we not have the coordinates?
  Wemissloc <- unique(Gen.and.loc[which(is.na(Gen.and.loc$Longitude)),][c("Code", "Year")])
  Wemissloc$We_miss <- "Coordinates missing"
  
  # For which code and year do we not have the HI?
  WemissHI <- unique(Gen.and.loc[which(is.na(Gen.and.loc$HI)),][c("Code", "Year")])
  WemissHI$We_miss <- "HI missing"
  
  # To write out :
  # write.csv(x = rbind(Wemissloc, WemissHI), file = "../output_data/information_missing.csv", row.names = FALSE)
  
  # Then correct :
  Gen.and.loc <- na.omit(Gen.and.loc)
  
  # Last check duplicates 
  which(duplicated(Gen.and.loc$Mouse_ID))
  
  # Print :
  Gen.and.loc
}

##
buildmap <- function(input, GenDF = GenDF, size = 2){
  
  # Merge with GenDF
  mergedDF <- merge(GenDF, input, by = "Mouse_ID")
  mergedDF$HI <- as.numeric(mergedDF$HI)
  
  # get a map
  margin <- 2
  
  area <- get_map(location =
                    c(min(mergedDF$Longitude - margin),
                      min(mergedDF$Latitude - margin),
                      max(mergedDF$Longitude + margin),
                      max(mergedDF$Latitude + margin)),
                  source = "stamen", maptype="toner-lite",
                  zoom = 7)
  
  #plot the map :
  ggmap(area) +
    geom_point(data = mergedDF, shape = 21, size = size,
               aes(Longitude, Latitude, fill = HI), alpha = 0.5) + # set up the points
    scale_fill_gradient("Hybrid\nindex", high="red",low="blue")   # set up the HI colors
}                  

## Function creating haplotypes:
myBioutifulHaplo <- function(input, GenDF){
  d <- ape::read.dna(input, format='fasta')
  
  ## Turn the row names into corresponding HI :
  seqnames <- data.frame(Mouse_ID = labels(d))
  
  matchname <- merge(seqnames, GenDF, by = "Mouse_ID", all.x = TRUE, sort = FALSE)
  
  matchname[which(is.na(matchname$HI)),]$HI <- "unknow_yet"
  
  match <- merge(seqnames, matchname, by = "Mouse_ID", sort = FALSE)
  
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
  mycols <- c(colorRampPalette(c("blue", "red"))(ncol(ind.hap) - 1), "grey")
  
  plot(net, size = attr(net, "freq"), pie = ind.hap, fast = TRUE, scale.ratio = 0.5, 
       cex = 1, bg = mycols)
  legend(-40, 40, colnames(ind.hap), fill = mycols, pch=19, ncol=2, cex = 0.6)
}
