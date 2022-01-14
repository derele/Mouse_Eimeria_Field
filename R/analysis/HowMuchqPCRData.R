library(tidyverse)

sota <- read.csv("data_products/SOTA_Data_Product.csv")

sota  %>% group_by(Year) %>%
    summarize(n = n(),              
              n.na = sum(is.na(delta_ct_cewe_MminusE)),
              n.Nna = sum(!is.na(delta_ct_cewe_MminusE)),
              n.MC = sum(!is.na(MC.Eimeria)),
              mD = mean(delta_ct_cewe_MminusE, na.rm=TRUE),
              mDINF = mean(delta_ct_cewe_MminusE[as.logical(MC.Eimeria)], na.rm=TRUE))


sota


sota  %>% group_by(Year) %>%
    summarize(n = n(),              
              OOn.na = sum(is.na(OPG)),
              OOn.Nna = sum(!is.na(OPG)),
              FEC.Nna = sum(!is.na(Feces_Weight)), 
              n.MC = sum(!is.na(MC.Eimeria)),
              mD = mean(delta_ct_cewe_MminusE, na.rm=TRUE),
              mDINF = mean(delta_ct_cewe_MminusE[as.logical(MC.Eimeria)], na.rm=TRUE))

sota
