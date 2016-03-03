convert.deg <-function(c){                                                                                                                                                                     
    c <- gsub(° , °,  c)                                                                                                                                                                   
    c <- gsub(' , .,  c)                                                                                                                                                                  
    x <- gsub(", ,  c)                                                                                                                                                                    
    z <- sapply((strsplit(x, [°\.])), as.numeric)                                                                                                                                           
    z[1, ] + z[2, ]/60 + z[3, ]/3600                                                                                                                                                           
} 
