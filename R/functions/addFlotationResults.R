addFlotationResults <- function(aDataFrame){
  ## Load data from oocysts counting 
  flotDF <- read.csv("../raw_data/Eimeria_detection/FINALOocysts2015to2017.csv")
  flotDF$OPG <- as.numeric(as.character(flotDF$OPG))
  flotDF <- flotDF[!is.na(flotDF$OPG),]
  
  ## Lorenzo count (in 1mL dilution) for comparison
  LorenzoDF <- read.csv("../raw_data/Eimeria_detection/Eimeria_oocysts_2015&2017_Lorenzo.csv")
  LorenzoDF <- LorenzoDF[!is.na(LorenzoDF$OPG),]
  
  ### Plot comparative Alice (dilution 0.1mL for most samples) and Lorenzo (dilution 1mL)
  compData <- merge(flotDF, LorenzoDF, by = "Mouse_ID", all = T)
  
  # merge with current data
  aDataFrame = merge(aDataFrame, flotDF, all = T)
  
  return(list(compData = compData,
              newDF = aDataFrame))
}

comparisonFlot <- function(aDataFrame){
  res <- addFlotationResults(aDataFrame)
  
  # How many samples new were detected by decreasing the dilution?
  N1 <- sum(
    res$compData$OPG.x > 0 & 
      res$compData$OPG.y == 0, na.rm = T)
  N1
  
  adjrsq <- summary(lm(formula = 
                         res$compData$OPG.x ~ 
                         res$compData$OPG.y))$adj.r.squared
  return(list(N1 = N1, adjrsq = adjrsq))
}

plotCompData <- function(aDataFrame){
  res = addFlotationResults(aDataFrame)
  plot1 <- ggplot(
    res$compData, aes(x = OPG.x, y = OPG.y)) +
    geom_point(alpha = .5, size = 4) +
    coord_equal(ratio=1) +
    xlab("counted in 0.1ml") + 
    ylab("counted in 1ml") + 
    geom_smooth(method = "lm", se = FALSE, col = "red") +
    geom_abline(intercept = 0, slope = 1, linetype = 3) +
    scale_y_log10() + 
    scale_x_log10() +
    theme_bw()
  return(plot1)
}
