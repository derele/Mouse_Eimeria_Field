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

is.geo.consistent <-  function (x, granularity=0.0001, consistent=TRUE) {
    select <- max(x$Latitude)-min(x$Latitude)< granularity &
        max(x$Longitude)-min(x$Longitude)< granularity
   if (consistent){return(x[select, ]) }
   else{return(x[!select,])}
}
