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
  
  ## Cleaned locations
  cleaned_loc <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/output_data/all_clean_localities.csv")
  cleaned_loc <- cleaned_loc[-1]
  
  ## Add dereje data (NB : No mice associated, just localities)
  ludo <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/raw_data/dureje_et_al._supplementary_table_s1.csv",
                   skip=1)
  ludo$Latitude <- convert.deg(ludo$Latitude)
  ludo$Longitude <- convert.deg(ludo$Longitude)
  ludo <- ludo[!ludo$Longitude==0, ]
  ludo <- ludo[!ludo$Latitude==0, ]
  ludo$Loc. <- paste("Loc",ludo$Loc., sep="_")
  names(ludo)[names(ludo)%in%"Loc."] <- "Code"
  # Convert :
  ludo$Latitude <- convert.deg(ludo$Latitude)
  ludo$Longitude <- convert.deg(ludo$Longitude)
  ludo$Year <- "review derele"
  # All loc until 2016 :
  All_loc <- rbind.match.columns(cleaned_loc, ludo)
  
  ### Genotypes
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
  
  # Add the data more complete :
  genotypes.2016.and.some.previous <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/raw_data/HIforEH_May2017.csv") 
  
  #######
  Gen.almost.tot <- data.frame(Mouse_ID = genotypes.2016.and.some.previous$PIN, Code = genotypes.2016.and.some.previous$Code,
                               Year = genotypes.2016.and.some.previous$Year, HI = genotypes.2016.and.some.previous$HI,
                               HI_NLoci = genotypes.2016.and.some.previous$HI_NLoci)
  Gen.almost.tot$Mouse_ID <- gsub(pattern = "SK", replacement = "SK_",x = Gen.almost.tot$Mouse_ID)
  
  Gen.tot <- merge(Gen.almost.tot, Gen14.15, by = c("Mouse_ID", "Code", "Year"), all = TRUE)
  
  # By default :
  Gen.tot$HI <- Gen.tot$HI.y
  Gen.tot$HI_NLoci <- "HI 6"
  
  # But if "better" HI calculated :
  Gen.tot[which(is.na(Gen.tot$HI)),]$HI_NLoci <- as.character(Gen.tot[which(is.na(Gen.tot$HI)),]$HI_NLoci.x)
  Gen.tot[which(is.na(Gen.tot$HI)),]$HI <- Gen.tot[which(is.na(Gen.tot$HI)),]$HI.x
  
  # Final :
  Gen.tot <- Gen.tot[c(1,2,3,8,9)]
  
  ## remove duplicates
  Gen.tot <- unique(Gen.tot)

  Gen.tot$HI <- round(as.numeric(as.character(Gen.tot$HI)), 2)
 
  # Add latitude / longitude
  Gen.tot <- merge(Gen.tot, All_loc, by = "Code", all.x = TRUE)
  
  # Print :
  Gen.tot
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