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
get.HI <- function (x){
  dom <- nchar(x) - nchar(gsub("d", "", x))
  mus <- nchar(x) - nchar(gsub("m", "", x))
  mus/(mus + dom)
}

## Calculate the hybrid index in a dataframe
get.HI.full <- function(df, markers.col = c("X332", "X347", "X65", "Tsx", "Btk",  "Syap1")){
  ## collapse m or d alleles in 1 column
  df$collapsed.GT <- apply(df[, markers.col],1, paste, collapse = "/")
  ## calculate hybrid index per mouse
  df$HI.calculated <- sapply(df$collapsed.GT, get.HI)
  return(df)
}

## Test if 2 localities form a cluster

## pairwise comparison of all localities 
pairwise.cluster.loc <- function (d, granularity=0.001){
    is.loc.cluster <- function(locA, locB){
        diff <- as.numeric(as.character(locA[c("Latitude", "Longitude")])) -
            as.numeric(as.character(locB[c("Latitude", "Longitude")]))
        ## both east-west and north-south are closer than granularity
        all (abs(diff) < granularity)
    }
    loc.combis <- combn(nrow(d), 2)
    long <- apply(loc.combis, 2, function(x)
        is.loc.cluster(d[x[1], ], d[x[2], ]))
    mat <- matrix(NA, nrow=nrow(d), ncol=nrow(d))
    mat[lower.tri(mat)] <- long
    mat <- t(mat)
    ## making the matrix symmetrical accross the diagonal
    diag(mat) <- 1
    mat[lower.tri(mat)] <- long
    mat
}


#### Problem: we need to recursivel cluster to cluster merge localites
#### which are not directly connected by via a thrid 
## cluster.loc <- function(d){
##     ## a list of lists to store each elements matches
##     clust.list <- apply(d, 1, function (x) which(x>0))
##     easy.clust.idx <- which(unlist(lapply(clust.list, length))==1)
##     to.clust.list<- clust.list[unlist(lapply(clust.list, length))>1]
##     combn(to.clust.list, 2, function (x) any(x[[1]]%in%x[[2]]), simplify=FALSE)
## }





# foo <- pairwise.cluster.loc(Trap)

## how many locs are sampled multiple times
# table(apply(foo, 2, sum))

## Convert distance in latitude/longitude to meters
measure <- function(lon1,lat1,lon2,lat2) {
  R <- 6378.137                                # radius of earth in Km
  dLat <- (lat2-lat1)*pi/180
  dLon <- (lon2-lon1)*pi/180
  a <- sin((dLat/2))^2 + cos(lat1*pi/180)*cos(lat2*pi/180)*(sin(dLon/2))^2
  c <- 2 * atan2(sqrt(a), sqrt(1-a))
  d <- R * c
  return (d * 1000)                            # distance in meters
}

## draw a map with colors depending on HI
HI.map <- function(df, size = 3, margin = 2, zoom = 7, alpha = 0.5){
  # get a map
  area <- get_map(location =
                    c(min(df$Longitude - margin),
                      min(df$Latitude - margin),
                      max(df$Longitude + margin),
                      max(df$Latitude + margin)),
                  source = "stamen", maptype="toner-lite",
                  zoom = zoom)
  #plot the map :
  ggmap(area) +
    geom_point(data = df, shape = 21, size = size,
               aes(Longitude, Latitude, fill = HI), alpha = alpha) + # set up the points
    scale_fill_gradient("Hybrid\nindex", high="red",low="blue")   # set up the HI colors
}                  

## Function creating haplotypes:
myBioutifulHaplo <- function(myfasta, dfwithHI){
  d <- ape::read.dna(myfasta, format='fasta')
  
  ## Turn the row names into corresponding HI :
  seqnames <- data.frame(Mouse_ID = labels(d))
  matchname <- merge(seqnames, dfwithHI, by = "Mouse_ID", all.x = TRUE, sort = FALSE)
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
