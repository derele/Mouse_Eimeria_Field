library(tidyverse)
library(visdat)
library(data.table)
library(stringr)



basics          <- c("Mouse_ID", "Sex", "Longitude", 
                     "Latitude", "Year", "HI", "HI_NLoci")

gen.loci        <- c("mtBamH", "YNPAR", "X332", "X347", "X65", "Tsx", "Btk", "Syap1",
                     "Es1", "Gpd1", "Idh1", "Mpi", "Np", "Sod1", "Es1C", "Gpd1C",
                     "Idh1C", "MpiC", "NpC", "Sod1C", "HI_NLoci",
                     "HI", "Zfy2", "Y", "Region")

dissection.cols <- c("Body_Weight", "Body_Length", "Tail_Length", "Status", "Spleen", 
                     "Left_Testis", "Right_Testis", "Seminal_Vesicles_Weight", "Liver",
                     "Sperm", "Left_Epididymis", "Right_Epididymis", 
                     "Right_Ovarium_Weight", "Left_Ovarium_Weight",
                     "Left_Embryo", "Right_Embryo", "Fleas", "Ticks", "Ectoparasites_Logical", 
                     "Arrival", "Dissection_Date", "Trap_Date", "Host")

Gene.Exp.cols   <- c("IFNy",        "CD4",         "Treg",   
                     "Treg17",      "Th1",         "Th17",   
                     "CD8",         "Act_CD8",     "IFNy_CD4",    "IL17A_CD4",  
                     "IFNy_CD8",    "IL.12",       "IRG6",        "CXCR3",           
                     "IL.6"    ,    "GBP2")


Crypto_qPCR.cols  <- c("Ct_mean", "Ct_mean_Ep", "Ct_mean_ABI", 
                       "Machine", "Measurements", "Tested_by", 
                       "qPCR_Date", "Oocyst_Predict", "Crypto_Positive")


Crypto_DNA.cols   <- c("ILWE_DNA_Content_ng.microliter", "ILWE_Tissue_used_up")



# add Data Product with filter for Samples that follow the 'AA_' pattern
    SOTA <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/SOTA_Data_Product.csv") %>% select(-X)
    SOTA <- SOTA[SOTA$Mouse_ID %like% "AA_", ]
    SOTA <- SOTA[colnames(SOTA) %in% c(basics)]
  
### add Crypto DNA Extraction Data
    Crypto_DNA   <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Cryptosporidium/DNA_Extraction_ILWE_2018_2021.csv")

### add Crypto_qPCR Data and filter out duplicates
    Crypto_qPCR           <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Cryptosporidium/Crypto_qPCR_2016_2019.csv") 


### calculate Oocysts with prediction model
    ABI_Best_thSC     <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Cryptosporidium/Crypto_Standard_Curve.csv")
      ABI_Best_thSC     <-  filter(ABI_Best_thSC, Ct_mean > 0)
      linear_model0     <- lm(log2(Amount_Oocysts) ~ Ct_mean, data = ABI_Best_thSC)
      Oocyst_Predict    <- 2^predict(linear_model0, newdata = Crypto_qPCR)
      Crypto_qPCR <- data.frame(Crypto_qPCR, Oocyst_Predict)
      Crypto_qPCR <- Crypto_qPCR %>% mutate(Oocyst_Predict = replace(Oocyst_Predict, Oocyst_Predict == "4292821751815.77", "0"))
      Crypto_qPCR$Oocyst_Predict <- as.integer(Crypto_qPCR$Oocyst_Predict)
    
## add Status (Crypto-positive or negative)
    Crypto_qPCR <- Crypto_qPCR %>% mutate(Crypto_Positive = ifelse(Ct_mean > 0, T, F))
  
### merging qPCR and and Extraction data
    Crypto_Detection  <- full_join(Crypto_qPCR[colnames(Crypto_qPCR) %in% c(basics, Crypto_qPCR.cols, Crypto_DNA.cols)], 
                         Crypto_DNA[colnames(Crypto_DNA) %in% c(basics, Crypto_qPCR.cols, Crypto_DNA.cols)], by = "Mouse_ID")

## add HI Data with SOTA, filter Crypto-positive Samples
    Crypto_Detection  <- full_join(Crypto_Detection[colnames(Crypto_Detection) %in% c("Mouse_ID", Crypto_qPCR.cols, Crypto_DNA.cols)], SOTA[colnames(SOTA) %in% c(basics)])


## add Columns
    # mus_caught
    # Crypto_mus_caught
    # Top_Location
    # Infection_Rate
    Crypto_pull <- Crypto_Detection %>% count(Longitude)%>% mutate(mus_caught = n) %>% select(Longitude, mus_caught) %>% arrange(mus_caught)
    Crypto_pull_pos <- Crypto_Detection %>% filter(Ct_mean > 0) %>% count(Longitude) %>% mutate(Crypto_mus_caught = n) %>% select(Longitude, Crypto_mus_caught) %>% arrange(Crypto_mus_caught)
    Crypto_Detection_1 <- left_join(Crypto_Detection, Crypto_pull)
    Crypto_Detection <- left_join(Crypto_Detection_1, Crypto_pull_pos)
    Crypto_Detection <- Crypto_Detection %>% replace_na(list(Crypto_mus_caught = 0))
    
    Crypto_Detection <- Crypto_Detection %>%
      mutate(Top_Location = Crypto_mus_caught >= 3,
             Infection_Rate = Crypto_mus_caught / mus_caught)
## write csv
    write.csv(Crypto_Detection, "data_products/Crypto_Detection.csv")
    