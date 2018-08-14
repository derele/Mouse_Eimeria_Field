addPCRresults <- function(aDataFrame, pathtodata = "../raw_data/Eimeria_detection/Inventory_contents_all.csv"){
  PCRdf <- read.csv(pathtodata)
  
  #correct wrong names
  toremove <- paste0(paste0("X",1:30, "_"), collapse = "|")
  names(PCRdf) <- gsub(toremove, "", names(PCRdf))
  
  names(PCRdf)[names(PCRdf)%in%"ID_mouse"] <- "Mouse_ID"
  names(PCRdf)[names(PCRdf)%in%"Year"] <- "year"
  
  PCRdf$Mouse_ID = gsub(" ", "", PCRdf$Mouse_ID) # fix the extra space
  
  # by default, I enter PCRstatus as negative, then overwrite
  PCRdf$PCRstatus = "negative"
  
  # PCR positive = one of the 3 other markers than AP5 sequenced 
  # (Ap5 was used for detection only, the other markers for confirmation)
  PCRdf$PCRstatus[PCRdf$`18S_Seq` == "positive" | 
                    PCRdf$COI_Seq == "positive" | 
                    PCRdf$ORF470_Seq == "positive"] <- "positive"
  
  # PCRstatus is NA if everything is NA
  PCRdf$PCRstatus[is.na(PCRdf$Ap5_PCR) & 
                    is.na(PCRdf$`18S_Seq`) &
                    is.na(PCRdf$COI_Seq) &
                    is.na(PCRdf$ORF470_Seq)] <- NA
  
  # merge with actual df
  aDataFrame <- merge(aDataFrame, PCRdf, by = "Mouse_ID", all = T)
  
  return(aDataFrame)
}