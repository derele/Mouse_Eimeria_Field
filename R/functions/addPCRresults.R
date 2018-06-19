addPCRresults <- function(aDataFrame){
  PCRdf <- read.csv("../raw_data/Eimeria_detection/Inventory_contents_all.csv")
  #correct wrong names
  toremove <- paste0(paste0("X",1:30, "_"), collapse = "|")
  names(PCRdf) <- gsub(toremove, "", names(PCRdf))
  
  names(PCRdf)[names(PCRdf)%in%"ID_mouse"] <- "Mouse_ID"
  
  # number of TRUE Ap5 positive (+ another marker sequenced successfully)
  PCRdf$oneSeqPositive <- PCRdf$`18S_Seq` == "positive" | 
    PCRdf$COI_Seq == "positive" | 
    PCRdf$ORF470_Seq == "positive"
  
  # positive with Ap5 and at least 1 sequence
  PCRdf$PCR.positive <-
    PCRdf$Ap5_PCR %in% "positive" & PCRdf$oneSeqPositive %in% TRUE
  
  # not tested yet
  PCRdf$PCR.positive[is.na(PCRdf$Ap5_PCR)] <- NA
  
  # merge with actual df
  aDataFrame <- merge(aDataFrame, PCRdf, all.x = T)
  
  return(aDataFrame)
}
