addqPCRresults <- function(aDataFrame, 
                           pathtoqPCR2016 = "../raw_data/Eimeria_detection/qPCR_2016.csv", 
                           pathtoqPCR2017 = "../raw_data/Eimeria_detection/qPCR_2017.csv"){

  ## Import data 2016
  qpcrData2016 <- read.csv(pathtoqPCR2016)
  names(qpcrData2016)[1] <- "Mouse_ID"
  
  ## Import data 2017
  qpcrData2017 <- read.csv(pathtoqPCR2017)
  
  ## merge all years
  qpcrData2016
  
  qpcrData2017Clean <- qpcrData2017["Mouse_ID"]
  
  # Add CEWE
  qpcrData2017Clean <- merge(qpcrData2017Clean, 
                             qpcrData2017[qpcrData2017$tissue == "CEWE", c("Mouse_ID", "deltaCt")],
                             all.x = T)
  names(qpcrData2017Clean)[names(qpcrData2017Clean) == "deltaCt"] <- "delta_ct_cewe"
  
  # Add ILWE
  qpcrData2017Clean <- merge(qpcrData2017Clean, 
                             qpcrData2017[qpcrData2017$tissue == "ILWE", c("Mouse_ID", "deltaCt")],
                             all.x = T)
  names(qpcrData2017Clean)[names(qpcrData2017Clean) == "deltaCt"] <- "delta_ct_ilwe"
  
  # Add observer
  qpcrData2017Clean$observer_qpcr <- "Lorenzo"
  
  # Merge both years
  qpcrData <- rbind(qpcrData2016, qpcrData2017Clean)
  
  #####
  
  # Did Enas calculate the other way around? 
  qpcrData$delta_ct_cewe[qpcrData$observer_qpcr == "Enas"] <- 
    - qpcrData$delta_ct_cewe[qpcrData$observer_qpcr == "Enas"]
  
  qpcrData$delta_ct_ilwe[qpcrData$observer_qpcr == "Enas"] <- 
    - qpcrData$delta_ct_ilwe[qpcrData$observer_qpcr == "Enas"]
  
  # deltaCT = ct eimeria - ct mouse. If high infection, low deltaCT
  # -deltaCT = ct mouse - ct eimeria
  qpcrData$qPCRsummary[qpcrData$delta_ct_cewe > 6 & qpcrData$delta_ct_ilwe > 6] <- "non infected"
  qpcrData$qPCRsummary[qpcrData$delta_ct_cewe < 6 & qpcrData$delta_ct_ilwe > 6] <- "infected cecum"
  qpcrData$qPCRsummary[qpcrData$delta_ct_cewe > 6 & qpcrData$delta_ct_ilwe < 6] <- "infected ileum"
  
  qpcrData$qPCRsummary[
    qpcrData$delta_ct_cewe < 6 & 
      qpcrData$delta_ct_ilwe < 6 & 
      qpcrData$delta_ct_cewe < qpcrData$delta_ct_ilwe] <- "cecum stronger"
  qpcrData$qPCRsummary[
    qpcrData$delta_ct_cewe < 6 & 
      qpcrData$delta_ct_ilwe < 6 & 
      qpcrData$delta_ct_cewe > qpcrData$delta_ct_ilwe] <- "ileum stronger"
  
  # Infected or not?
  qpcrData$qPCRstatus <- "positive"
  qpcrData$qPCRstatus[is.na(qpcrData$qPCRsummary)] <- NA
  qpcrData$qPCRstatus[qpcrData$qPCRsummary %in% "non infected"] <- "negative"
  
  # and keep the infected segment value OR the higher value 
  qpcrData$delta_ct[
    qpcrData$qPCRsummary %in% c("infected cecum", "cecum stronger")] <- 
    qpcrData$delta_ct_cewe[
      qpcrData$qPCRsummary %in% c("infected cecum", "cecum stronger")] 
  
  qpcrData$delta_ct[
    qpcrData$qPCRsummary %in% c("infected ileum", "ileum stronger")] <- 
    qpcrData$delta_ct_ilwe[
      qpcrData$qPCRsummary %in% c("infected ileum", "ileum stronger")] 
  
  # Turn around
  qpcrData$delta_ct_MminusE <- - qpcrData$delta_ct
  
  # Set floor values
  qpcrData$delta_ct_MminusE[is.na(qpcrData$delta_ct_MminusE)] <- -6
  
  # To pass positive I add 6 to all
  # qpcrData$delta_ct_MminusE <- qpcrData$delta_ct_MminusE + 6

  # merge
  aDataFrame <- merge(aDataFrame, qpcrData, by = "Mouse_ID", all = T)
  
  return(aDataFrame)
}
