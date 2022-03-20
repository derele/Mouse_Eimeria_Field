library(tidyverse)


df_results <- read.csv('https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/qPCR_faeces_field_eimeria_22')

#### PIVOT AND UNIFICATION WITH OTHER HZ DATA ##################################

# correct colnames
colnames(df_results)[colnames(df_results)%in%"Sample"]  <- "Mouse_ID"
colnames(df_results)[colnames(df_results)%in%"Cq Mean"] <- "FEC_Eim_Ct" # specific to this assay, to avoid confusion with other qPCRs
# correct the HZ Mouse_IDs --> pattern AA_0XXX
df_results$Mouse_ID <- gsub(pattern = "A", replacement = "AA_0", x = df_results$Mouse_ID)
df_results$Mouse_ID <- gsub(pattern = "AA_0AA_", replacement = "AA_", x = df_results$Mouse_ID)
# correct the Standard Curve (SC) Mouse_IDs --> pattern VXX_X
df_results$Mouse_ID <- gsub(pattern = "v10_5", replacement = "V10_5", x = df_results$Mouse_ID)
df_results$Mouse_ID <- gsub(pattern = "v10_2", replacement = "V10_2", x = df_results$Mouse_ID)
# correct type of some cols to numeric 
df_results[,c(1, 10, 13:16, 18, 20:25)] <- lapply(df_results[, c(1, 10, 13:16, 18, 20:25)], as.numeric)
# correct task type (wild = UNKNOWN, SC = STANDARD, NTC = NTC)
df_results <- df_results %>% mutate(Task = case_when(Mouse_ID == 'NTC' ~ 'NTC',
                                                     Mouse_ID == 'NTCE' ~ 'NTC',
                                                     Mouse_ID == 'NTCM' ~ 'NTC',
                                                     Mouse_ID == 'V10_0' ~ 'STANDARD',
                                                     Mouse_ID == 'V10_1' ~ 'STANDARD',
                                                     Mouse_ID == 'V10_2' ~ 'STANDARD',
                                                     Mouse_ID == 'V10_3' ~ 'STANDARD',
                                                     Mouse_ID == 'V10_4' ~ 'STANDARD',
                                                     Mouse_ID == 'V10_5' ~ 'STANDARD',
                                                     Mouse_ID == 'V10_6' ~ 'STANDARD',
                                                     is.na(Task) ~ 'UNKNOWN'))
# specify Ct measurements from well plate position
df_results <- df_results %>% mutate(Ct_index = case_when (Well.Position %in% c('A1', 'A4', 'A7', 'A10', 
                                                                                 'B1', 'B4', 'B7', 'B10', 
                                                                                 'C1', 'C4', 'C7', 'C10', 
                                                                                 'D1', 'D4', 'D7', 'D10',
                                                                                 'E1', 'E4', 'E7', 'E10',
                                                                                 'F1', 'F4', 'F7', 'F10',
                                                                                 'G1', 'G4', 'G7', 'G10',
                                                                                 'H1', 'H4', 'H7', 'H10') ~ 'Ct1', 
                                                          Well.Position %in% c('A2', 'A5', 'A8', 'A11', 
                                                                                 'B2', 'B5', 'B8', 'B11', 
                                                                                 'C2', 'C5', 'C8', 'C11', 
                                                                                 'D2', 'D5', 'D8', 'D11',
                                                                                 'E2', 'E5', 'E8', 'E11',
                                                                                 'F2', 'F5', 'F8', 'F11',
                                                                                 'G2', 'G5', 'G8', 'G11',
                                                                                 'H2', 'H5', 'H8', 'H11')  ~ 'Ct2', 
                                                          Well.Position %in% c('A3', 'A6', 'A9', 'A12', 
                                                                                 'B3', 'B6', 'B9', 'B12', 
                                                                                 'C3', 'C6', 'C9', 'C12', 
                                                                                 'D3', 'D6', 'D9', 'D12',
                                                                                 'E3', 'E6', 'E9', 'E12',
                                                                                 'F3', 'F6', 'F9', 'F12',
                                                                                 'G3', 'G6', 'G9', 'G12',
                                                                                 'H3', 'H6', 'H9', 'H12')  ~ 'Ct3')) 


# pivot table to have all measurements per Mouse_ID in one row (not three)
df_results <- pivot_wider(df_results, names_from = "Ct_index", values_from = "Cq") %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% 
  fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T)
# make new cols numeric
df_results$FEC_Eim_Ct1 <- as.numeric(df_results$Ct1);
df_results$FEC_Eim_Ct2 <- as.numeric(df_results$Ct2);
df_results$FEC_Eim_Ct3 <- as.numeric(df_results$Ct3)

# (Un)Select Columns for final Data Product:
df_results <- df_results %>% select(-c(Well, `Well.Position`, Omit, Quencher, `Curve.Quality`, `Result.Quality.Issues`, `Auto.Threshold`))


#write the data frame in a csv file 
write.csv(df_results, "~/Documents/Mouse_Eimeria_Field/data_products/FEC_EqPCR_DataProduct.csv", row.names=FALSE)


