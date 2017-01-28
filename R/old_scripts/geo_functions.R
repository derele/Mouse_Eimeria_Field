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


get.HIX <- function (gt.data.cols, rows=FALSE){
    collapsed.GT <- apply(gt.data.cols, 1, paste, collapse = "/")
    if (rows){ ## to merge over mice from e.g. one locality (different rows)
        collapsed.GT <- paste(collapsed.GT, collapse = "/")
    }
    dom <- nchar(collapsed.GT) - nchar(gsub("d", "", collapsed.GT))
    mus <- nchar(collapsed.GT) - nchar(gsub("m", "", collapsed.GT))
    mus/(mus + dom)
}
}
