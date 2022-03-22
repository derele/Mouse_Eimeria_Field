library(tidyverse)


# STATUS: There are a lot of samples without a sample name.. I assume plate setup error led to this.
# Please double-check the data in your workbooks to ensure that at least some of the nameless samples that were positive in qPCR
# can be assigned to the according Mouse_ID (not many samples, but still - we have the data but cannot connect it..)



CEWE_Clean <- read.csv('https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/qPCR_CEWE_eimeria_21_22.csv')

# correct colnames
colnames(CEWE_Clean)[colnames(CEWE_Clean)%in%"Sample"]  <- "Mouse_ID"
colnames(CEWE_Clean)[colnames(CEWE_Clean)%in%"Cq"]      <- "Ct"
colnames(CEWE_Clean)[colnames(CEWE_Clean)%in%"Cq.Mean"] <- "Ct_Mean"
colnames(CEWE_Clean)[colnames(CEWE_Clean)%in%"Cq.SD"]   <- "Ct_StDev"

# correct the HZ Mouse_IDs --> pattern AA_0XXX
CEWE_Clean$Mouse_ID <- gsub(pattern = "A", replacement = "AA_0", x = CEWE_Clean$Mouse_ID)
CEWE_Clean$Mouse_ID <- gsub(pattern = "AA_0AA_", replacement = "AA_", x = CEWE_Clean$Mouse_ID)

# correct type of some cols to numeric 
# currently leads to NA for undetermined Ct values --> maybe 0 is better to add information that the sample was run in qPCR 
# (Oocyst prediction (if necessary) could be corrected)
CEWE_Clean[,c(1, 10, 13:16, 18, 20:25)] <- lapply(CEWE_Clean[, c(1, 10, 13:16, 18, 20:25)], as.numeric)

# correct Target spelling differences
CEWE_Clean$Target[grep("eim*.", CEWE_Clean$Target)] <- "Eim"
CEWE_Clean$Target[grep("mus*.", CEWE_Clean$Target)] <- "Mus"

# remove (unnecessary) columns --> please correct me in case some are important for YOU!
#CEWE_Clean <- CEWE_Clean %>% select(-c(Well, Omit, Reporter, Quencher, Amp.Score, Amp.Status, Curve.Quality, Result.Quality.Issues, Cq.Confidence, Auto.Threshold, Threshold, Auto.Baseline, Baseline.Start, Baseline.End))

# subset all samples that target Mus
subM <- CEWE_Clean$Mouse_ID[grep('*M', CEWE_Clean$Mouse_ID)]
subM <- CEWE_Clean %>% filter(Mouse_ID %in% subM)
colnames(subM) <- paste(colnames(subM), 'Mus', sep = '_')
colnames(subM)[colnames(subM) %in% c( "Well.Position_Mus", "Mouse_ID_Mus", 'Target_Mus', "Task_Mus", "Ct_Mean_Mus", 'plate_Mus')]       <- c("Well.Position", "Mouse_ID", 'Target', "Task", "Ct_Mean", 'plate')
# add Ct_index for pivot
subM <- subM %>% mutate(Ct_index = case_when (Well.Position %in% c('A4', 'B4', 'C4', 'D4', 'E4', 'F4', 'G4', 'H4',
                                                                   'A10', 'B10', 'C10', 'D10', 'E10', 'F10', 'G10', 'H10') ~ 'Ct1_Mus',
                                              Well.Position %in% c('A5', 'B5', 'C5', 'D5', 'E5', 'F5', 'G5', 'H5',
                                                                   'A11', 'B11', 'C11', 'D11', 'E11', 'F11', 'G11', 'H11') ~ 'Ct2_Mus',
                                              Well.Position %in% c('A6', 'B6', 'C6', 'D6', 'E6', 'F6', 'G6', 'H6',
                                                                   'A12', 'B12', 'C12', 'D12', 'E12', 'F12', 'G12', 'H12') ~ 'Ct3_Mus'),
                        Target = case_when(Target == 'Mus'   ~ 'Mus',
                                           Target == 'Eim'   ~ 'Mus',
                                           Target == 'CDC42' ~ 'CDC42',
                                           Target == 'COI'   ~ 'COI'))

# subset all samples that target Eim
subE <- CEWE_Clean$Mouse_ID[grep('*E', CEWE_Clean$Mouse_ID)]
subE <- CEWE_Clean %>% filter(Mouse_ID %in% subE)
colnames(subE) <- paste(colnames(subE), 'Eim', sep = '_')
colnames(subE)[colnames(subE) %in% c( "Well.Position_Eim", "Mouse_ID_Eim", 'Target_Eim', "Task_Eim", "Ct_Mean_Eim", 'plate_Eim')]       <- c("Well.Position", "Mouse_ID", 'Target', "Task", "Ct_Mean", 'plate')
# add Ct_index for pivot
subE <- subE %>% mutate(Ct_index = case_when (Well.Position %in% c('A1', 'B1', 'C1', 'D1', 'E1', 'F1', 'G1', 'H1',
                                                                   'A7', 'B7', 'C7', 'D7', 'E7', 'F7', 'G7', 'H7') ~ 'Ct1_Eim',
                                              Well.Position %in% c('A2', 'B2', 'C2', 'D2', 'E2', 'F2', 'G2', 'H2',
                                                                   'A8', 'B8', 'C8', 'D8', 'E8', 'F8', 'G8', 'H8') ~ 'Ct2_Eim',
                                              Well.Position %in% c('A3', 'B3', 'C3', 'D3', 'E3', 'F3', 'G3', 'H3',
                                                                   'A9', 'B9', 'C9', 'D9', 'E9', 'F9', 'G9', 'H9') ~ 'Ct3_Eim'),
                        Target = case_when(Target == 'Mus'   ~ 'Eim',
                                           Target == 'Eim'   ~ 'Eim',
                                           Target == 'COI'   ~ 'COI'))


# merge the subsets
CEWE_Clean <- full_join(subM, subE) %>% arrange(Mouse_ID)

# correct Mouse_ID (didn't do that earlier because I used it for subsetting)
CEWE_Clean$Mouse_ID <- gsub(pattern = "*M", replacement = "", x = CEWE_Clean$Mouse_ID)
CEWE_Clean$Mouse_ID <- gsub(pattern = "*E", replacement = "", x = CEWE_Clean$Mouse_ID)



# pivot table to have all measurements per Mouse_ID in one row (not three)
CEWE_Clean <- pivot_wider(CEWE_Clean, names_from = "Ct_index", values_from = 'Ct_Mean') %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% 
  fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) %>% filter(Task != 'NTC')

# (Un)Select Columns for final Data Product:
CEWE_Clean <- CEWE_Clean %>% 
  mutate(delta_ct_cewe_MminusE = (Ct_Mus - Ct_Eim)) %>%
  select(Mouse_ID, delta_ct_cewe_MminusE, 
         Ct_Mus, Ct1_Mus, Ct2_Mus, Ct3_Mus, 
         Ct_Eim, Ct1_Eim, Ct2_Eim, Ct3_Eim, 
         Tm1_Mus, Tm2_Mus, Tm3_Mus, Tm4_Mus, 
         Tm1_Eim, Tm2_Eim, Tm3_Eim, Tm4_Eim, 
         Ct_StDev_Mus, Ct_StDev_Eim, plate)

#write the data frame in a csv file 
write.csv(CEWE_Clean, "~/Documents/Mouse_Eimeria_Field/data_products/HZ21_CEWE_EqPCR_DataProduct.csv", row.names=FALSE)


