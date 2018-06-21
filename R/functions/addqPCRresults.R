addqPCRresults <- function(aDataFrame){
  ## Import data 2016
  qpcrData <- read.csv("../raw_data/Eimeria_detection/qPCR_2016.csv")
  names(qpcrData)[1] <- "Mouse_ID"
  
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
