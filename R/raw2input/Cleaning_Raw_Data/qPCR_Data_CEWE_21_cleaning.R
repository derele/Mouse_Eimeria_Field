library(xlsx)
library(tidyr)
library(dplyr)
library(janitor)
library(readr)

#creating script to read all the results file produced by the programm QuantStudio 1
#the result files are located at: GitHub/Mouse_Eimeria_Field/data_input/qPCR_CEWE_21/Results_Files/


#1 change the current working directory to the location where the qPCR data files are
setwd("~/GitHub/Mouse_Eimeria_Field/data_input/qPCR_CEWE_21/Results_Files/")

#2. Create a list of names, out of the files in this repository
#list those files
list_CEWE <- as.list(list.files()) #now we have a list of the data frame names

#We need to create a list out of the names 
list_names <- as.vector(unlist(list_CEWE))

#write a function to specify how to read the qPCR files
read_qPCR_file <- function(x) {
  
  df1 <- xlsx::read.xlsx(x, sheetIndex = 1)
  #get the file name of the file
  #the name of the file is the second column of this file
  #to get that name we can start by selecting this column
  #and then getting the name out of it
  filename <- colnames(df1[2])
  #remove unecessary rows of the data frame
  #everything before actual data
  df1 <- df1 %>%
    filter(!row_number() %in% c(1:23))
  #change the column names to the names of the first row
  colnames(df1) <- df1[1, ]
  #Now remove the first row
  df1 <- df1 %>% filter(!row_number() %in% 1)
  #add a column with the name of the plate 
  df1 <- df1 %>% mutate(plate = filename)
  
}

#apply the function you created in the last step to each of the elements (names of files)
#of the list list_names
#in this way you are creating a list of data frames from the result files produced by the 
#program of the machine
list_results <- lapply(list_names, read_qPCR_file)

#show the data frame consisting of each result data file
df_results <- Reduce(rbind, list_results)

#remove duplicates
df <- unique(df_results)

setwd("~/GitHub/Mouse_Eimeria_Field/")

# correct colnames
colnames(df)[colnames(df)%in%"Sample"]  <- "Mouse_ID"
colnames(df)[colnames(df)%in%"Cq"]      <- "Ct"
colnames(df)[colnames(df)%in%"Cq.Mean"] <- "Ct_Mean"
colnames(df)[colnames(df)%in%"Cq.SD"]   <- "Ct_StDev"

# correct the HZ Mouse_IDs --> pattern AA_0XXX
df$Mouse_ID <- gsub(pattern = "A", replacement = "AA_0", x = df$Mouse_ID)
df$Mouse_ID <- gsub(pattern = "AA_0AA_", replacement = "AA_", x = df$Mouse_ID)

# correct type of some cols to numeric 
# currently leads to NA for undetermined Ct values --> maybe 0 is better to add information that the sample was run in qPCR 
# (Oocyst prediction (if necessary) could be corrected)
df[,c(1, 10, 13:16, 18, 20:25)] <- lapply(df[, c(1, 10, 13:16, 18, 20:25)], as.numeric)

# correct Target spelling differences
df$Target[grep("eim*.", df$Target)] <- "Eim"
df$Target[grep("mus*.", df$Target)] <- "Mus"
df_preClean <- df

# Start Cleaning
CEWE_Clean <- df_preClean

# all samples that target Mus
CEWE_Clean_M <- CEWE_Clean$Mouse_ID[grep('*M', CEWE_Clean$Mouse_ID)]
CEWE_Clean_M <- CEWE_Clean %>% filter(Mouse_ID %in% CEWE_Clean_M)
colnames(CEWE_Clean_M) <- paste(colnames(CEWE_Clean_M), 'Mus', sep = '_')

CEWE_Clean_M <- dplyr::rename(CEWE_Clean_M, 
                              "Well.Position" = "Well Position_Mus",
                              "Mouse_ID" = "Mouse_ID_Mus",
                              "Target" = "Target_Mus", 
                              "Task" = "Task_Mus", 
                              "plate" = "plate_Mus")



# add Ct_index for pivot
CEWE_Clean_M <- CEWE_Clean_M %>% dplyr::mutate(Ct_index = case_when (Well.Position %in% c("A4", "B4", "C4", "D4", "E4", "F4", "G4", "H4",
                                                                                   "A10", "B10", "C10", "D10", "E10", "F10", "G10", "H10") ~ "Ct1_Mus",
                                                              Well.Position %in% c("A5", "B5", "C5", "D5", "E5", "F5", "G5", "H5",
                                                                                   "A11", "B11", "C11", "D11", "E11", "F11", "G11", "H11") ~ "Ct2_Mus",
                                                              Well.Position %in% c("A6", "B6", "C6", "D6", "E6", "F6", "G6", "H6",
                                                                                   "A12", "B12", "C12", "D12", "E12", "F12", "G12", "H12") ~ "Ct3_Mus",
                                                              TRUE ~ ""),
                                        Target = case_when(Target == "Mus"   ~ "Mus",
                                                           Target == "Eim"   ~ "Mus",
                                                           Target == "CDC42" ~ "CDC42",
                                                           Target == "COI"   ~ "COI",
                                                           TRUE ~ ""))


CEWE_Clean_M <- pivot_wider(CEWE_Clean_M, names_from = "Ct_index", values_from = "Ct_Mus") %>% 
  arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% 
  fill(c(everything()), .direction = "downup") %>% 
  ungroup() %>% 
  distinct(Mouse_ID, .keep_all = T) %>% 
  dplyr::filter(Task != "NTC")


# all samples that target Eim
CEWE_Clean_E <- CEWE_Clean$Mouse_ID[grep('*E', CEWE_Clean$Mouse_ID)]
CEWE_Clean_E <- CEWE_Clean %>% filter(Mouse_ID %in% CEWE_Clean_E)
colnames(CEWE_Clean_E) <- paste(colnames(CEWE_Clean_E), 'Eim', sep = '_')


CEWE_Clean_E <- dplyr::rename(CEWE_Clean_E, 
              "Well.Position" = "Well Position_Eim",
              "Mouse_ID" = "Mouse_ID_Eim",
              "Target" = "Target_Eim", 
              "Task" = "Task_Eim", 
              "plate" = "plate_Eim")



# add Ct_index for pivot
CEWE_Clean_E <- CEWE_Clean_E %>% 
  dplyr::mutate(Ct_index = case_when (Well.Position %in% c('A1', 'B1', 'C1', 'D1', 
                                                           'E1', 'F1', 'G1', 'H1',
                                                           'A7', 'B7', 'C7', 'D7', 
                                                           'E7', 'F7', 'G7', 'H7') ~ 'Ct1_Eim',
                                      Well.Position %in% c('A2', 'B2', 'C2', 'D2', 
                                                           'E2', 'F2', 'G2', 'H2',
                                                           'A8', 'B8', 'C8', 'D8', 
                                                           'E8', 'F8', 'G8', 'H8') ~ 'Ct2_Eim',
                                      Well.Position %in% c('A3', 'B3', 'C3', 'D3', 
                                                           'E3', 'F3', 'G3', 'H3',
                                                           'A9', 'B9', 'C9', 'D9', 
                                                           'E9', 'F9', 'G9', 'H9') ~ 'Ct3_Eim'),
                                        Target = case_when(Target == 'Mus'   ~ 'Eim',
                                                           Target == 'Eim'   ~ 'Eim',
                                                           Target == 'COI'   ~ 'COI'))

CEWE_Clean_E <- pivot_wider(CEWE_Clean_E, names_from = "Ct_index", values_from = 'Ct_Eim') %>% 
  arrange(Mouse_ID) %>% 
  group_by(Mouse_ID) %>% 
  fill(c(everything()), .direction = "downup") %>% 
  ungroup() %>% 
  distinct(Mouse_ID, .keep_all = T) %>% 
  filter(Task != 'NTC')


# merge the E and M sets
CEWE_Clean <- full_join(CEWE_Clean_M, CEWE_Clean_E)

# remove Mouse_ID add ons
CEWE_Clean$Mouse_ID <- gsub(pattern = "*M", replacement = "", x = CEWE_Clean$Mouse_ID)
CEWE_Clean$Mouse_ID <- gsub(pattern = "*E", replacement = "", x = CEWE_Clean$Mouse_ID)

# find, fill, remove duplicates
CEWE_Clean <- CEWE_Clean %>% 
  arrange(Mouse_ID) %>% 
  group_by(Mouse_ID) %>% 
  fill(c(everything()), .direction = "downup") %>% 
  ungroup() %>% 
  distinct(Mouse_ID, .keep_all = T) %>% 
  filter(Task != 'NTC')
rm(CEWE_Clean_M)
rm(CEWE_Clean_E)

# (Un)Select Columns for final Data Product:
colnames(CEWE_Clean)
CEWE_Clean <- dplyr::rename(CEWE_Clean, 
                            Ct_Mean_Mus = "Cq Mean_Mus",
                            Ct_Mean_Eim = "Cq Mean_Eim",
                            Ct_StDev_Mus = "Cq SD_Mus",
                            Ct_StDev_Eim = "Cq SD_Eim")


CEWE_Clean <- CEWE_Clean %>% 
  dplyr::mutate(delta_ct_cewe_MminusE = (Ct_Mean_Mus - Ct_Mean_Eim)) %>%
  dplyr::select(Mouse_ID, delta_ct_cewe_MminusE, 
         Ct_Mean_Mus, Ct1_Mus, Ct2_Mus, Ct3_Mus, 
         Ct_Mean_Eim, Ct1_Eim, Ct2_Eim, Ct3_Eim,
         Tm1_Mus, Tm2_Mus, Tm3_Mus, Tm4_Mus, 
         Tm1_Eim, Tm2_Eim, Tm3_Eim, Tm4_Eim, 
         Ct_StDev_Mus, Ct_StDev_Eim, plate)

# write WILD SAMPLES data frame in a csv file
write.csv(CEWE_Clean, "data_input/CEWE_EqPCR_Data_Product.csv", row.names=FALSE)

# STANDARD CURVE
SC <- df_preClean %>% filter(Task %in% c('STANDARD'))
SC_1 <- SC %>% filter(plate == 'qPCR_eimeria_field_CEWE_JK_04032022.eds')   %>% mutate(SC = 1)
SC_2 <- SC %>% filter(plate == 'qPCR_eimeria_field_CEWE_JK_04032022_2.eds') %>% mutate(SC = 2)
SC <- full_join(SC_1, SC_2); rm(SC_1); rm(SC_2)
colnames(SC) <- paste(colnames(SC), 'Eim', sep = '_')

SC <- dplyr::rename(SC, "Well.Position" = "Well Position_Eim", 
                    "Mouse_ID" = "Mouse_ID_Eim", 
                    "Target" =  "Target_Eim",
                    "Task" = "Task_Eim", 
                    "plate" = "plate_Eim", 
                    "SC" = "SC_Eim",
                    Ct_Mean_Eim = "Cq Mean_Eim",
                    Ct_StDev_Eim = "Cq SD_Eim")


SC$Target <- gsub(pattern = "CDC42", replacement = "COI", x = SC$Target)
SC <- SC %>% dplyr::mutate(Ct_index = case_when (Well.Position %in% c('A1', 'A4', 'A7', 'A10', 'B1', 'B4', 'B7', 'B10', 'C1', 'C4', 'C7', 'C10', 'D1', 'D4', 'D7', 'D10', 'E1', 'E4', 'E7', 'E10', 'F1', 'F4', 'F7', 'F10', 'G1', 'G4', 'G7', 'G10', 'H1', 'H4', 'H7', 'H10')  ~ 'Ct1_Eim', 
                                          Well.Position %in% c('A2', 'A5', 'A8', 'A11', 'B2', 'B5', 'B8', 'B11', 'C2', 'C5', 'C8', 'C11', 'D2', 'D5', 'D8', 'D11', 'E2', 'E5', 'E8', 'E11', 'F2', 'F5', 'F8', 'F11', 'G2', 'G5', 'G8', 'G11', 'H2', 'H5', 'H8', 'H11')  ~ 'Ct2_Eim', 
                                          Well.Position %in% c('A3', 'A6', 'A9', 'A12', 'B3', 'B6', 'B9', 'B12', 'C3', 'C6', 'C9', 'C12', 'D3', 'D6', 'D9', 'D12', 'E3', 'E6', 'E9', 'E12', 'F3', 'F6', 'F9', 'F12', 'G3', 'G6', 'G9', 'G12', 'H3', 'H6', 'H9', 'H12')  ~ 'Ct3_Eim')) 

SC <- pivot_wider(SC, names_from = "Ct_index", values_from = 'Ct_Eim') %>% arrange(Mouse_ID, SC) %>% group_by(Mouse_ID, SC) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, SC, .keep_all = T)

SC <- SC %>% dplyr::select(Mouse_ID, SC, Ct_Mean_Eim, Ct1_Eim, Ct2_Eim, Ct3_Eim, Tm1_Eim, Tm2_Eim, Tm3_Eim, Tm4_Eim, Ct_StDev_Eim, plate)

# write STANDARD CUFVE data frame in a csv file
write.csv(SC, "data_input/HZ21_CEWE_EqPCR_SC.csv", row.names=FALSE)
