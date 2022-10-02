library(xlsx)
library(tidyr)
library(dplyr)
library(janitor)
library(readr)

################################################################################
####  Eimeria qPCR  ############################################################
################################################################################


#create a list out of the names of each excel file (maybe you have to do that locally and upload it to github)
#eimeria_plates <- list.files(path = "~/Documents/Github/Mouse_Eimeria_Field/data_input/qPCR_FEC_22/results/eimeria_plates/")
#write.csv(eimeria_plates, "~/Documents/Github/Mouse_Eimeria_Field/data_input/qPCR_FEC_22/results/eimeria_plates.csv")
eimeria_plates      <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/qPCR_FEC_22/results/eimeria_plates.csv")
eimeria_plates_josi <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/qPCR_faeces_21/Results/Filenames_qpcr_results.csv")

#read the table whith the names of each file, call it NT (= Name table)
eim_list      <- as.list(eimeria_plates$x)
eim_list_josi <- as.list(eimeria_plates_josi$qPCR.Results.file.names)


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
  #different in Crypto (22) vs. Eimeria (23)!
  df1 <- df1 %>%
    filter(!row_number() %in% c(1:23))
  #change the column names to the names of the first row
  colnames(df1) <- df1[1, ]
  #Now remove the first row
  df1 <- df1 %>% filter(!row_number() %in% 1)
  #add a column with the name of the plate 
  df1 <- df1 %>% mutate(plate = filename)
  
}

setwd("/Users/finnlo/Documents/Github/Mouse_Eimeria_Field/data_input/qPCR_FEC_22/results/eimeria_plates/")
results <- lapply(eim_list, read_qPCR_file)

setwd("~/Documents//GitHub/Mouse_Eimeria_Field/data_input/qPCR_faeces_21/Results/Results_files/")
results_josi <- lapply(eim_list_josi, read_qPCR_file)

#show the data frame consisting of each result data file
eim_results      <- Reduce(rbind, results)
eim_results_josi <- Reduce(rbind, results_josi)

# combine Josi's data (HZ18-21) and HZ22
eim_results      <- bind_rows(eim_results, eim_results_josi)

#remove duplicates (eliminates some NTCs)
eim_results <- unique(eim_results)

# correct colnames
colnames(eim_results)[colnames(eim_results)%in%"Sample"]  <- "Mouse_ID"
colnames(eim_results)[colnames(eim_results)%in%"Cq Mean"] <- "FEC_Eim_Ct" # specific to this assay, to avoid confusion with other qPCRs

# correct the HZ Mouse_IDs --> pattern AA_0XXX
eim_results$Mouse_ID <- gsub(pattern = "AA", replacement = "AA_0", x = eim_results$Mouse_ID)
eim_results$Mouse_ID <- gsub(pattern = "AA_0_", replacement = "AA_", x = eim_results$Mouse_ID)

# correct the Standard Curve (SC) Mouse_IDs --> pattern VXX_X
eim_results$Mouse_ID <- gsub(pattern = "v10_5", replacement = "V10_5", x = eim_results$Mouse_ID)
eim_results$Mouse_ID <- gsub(pattern = "v10_2", replacement = "V10_2", x = eim_results$Mouse_ID)


# correct colname for easier selecting
eim_results <- dplyr::rename(eim_results, Well.Position = "Well Position",
                             Curve.Quality = "Curve Quality",
                             Result.Quality.Issues = "Result Quality Issues",
                             Auto.Threshold = "Auto Threshold",
                             Amp.Status = `Amp Status`,
                             Amp.Score = `Amp Score`,
                             Cq.Confidence = `Cq Confidence`,
                             Cq.SD = `Cq SD`)

# correct type of some cols to numeric 
eim_results[,c(1, 10, 13:16, 18, 20:25)] <- lapply(eim_results[, c(1, 10, 13:16, 18, 20:25)], as.numeric)

# correct row input of Cq --> undetermined to 0 
# a 0 means that the sample was run, but nothing was measured, 
# whereas NA would mean the absence of the sample..
# BEWARE: For oocyst prediction this 0 will lead to a super high prediction, must be corrected to 0 again!
eim_results <- eim_results %>% mutate(FEC_Eim_Ct = case_when(is.na(FEC_Eim_Ct) ~ 0,
                                                             !is.na(FEC_Eim_Ct) ~ FEC_Eim_Ct))
eim_results$FEC_Eim_Ct <- round(eim_results$FEC_Eim_Ct, 2)


# specify Ct measurements from well plate position
eim_results <- eim_results %>% mutate(Ct_index = case_when (Well.Position %in% c('A1', 'A4', 'A7', 'A10', 
                                                                                 'B1', 'B4', 'B7', 'B10', 
                                                                                 'C1', 'C4', 'C7', 'C10', 
                                                                                 'D1', 'D4', 'D7', 'D10',
                                                                                 'E1', 'E4', 'E7', 'E10',
                                                                                 'F1', 'F4', 'F7', 'F10',
                                                                                 'G1', 'G4', 'G7', 'G10',
                                                                                 'H1', 'H4', 'H7', 'H10') ~ 'FEC_Eim_Ct1', 
                                                            Well.Position %in% c('A2', 'A5', 'A8', 'A11', 
                                                                                 'B2', 'B5', 'B8', 'B11', 
                                                                                 'C2', 'C5', 'C8', 'C11', 
                                                                                 'D2', 'D5', 'D8', 'D11',
                                                                                 'E2', 'E5', 'E8', 'E11',
                                                                                 'F2', 'F5', 'F8', 'F11',
                                                                                 'G2', 'G5', 'G8', 'G11',
                                                                                 'H2', 'H5', 'H8', 'H11')  ~ 'FEC_Eim_Ct2', 
                                                            Well.Position %in% c('A3', 'A6', 'A9', 'A12', 
                                                                                 'B3', 'B6', 'B9', 'B12', 
                                                                                 'C3', 'C6', 'C9', 'C12', 
                                                                                 'D3', 'D6', 'D9', 'D12',
                                                                                 'E3', 'E6', 'E9', 'E12',
                                                                                 'F3', 'F6', 'F9', 'F12',
                                                                                 'G3', 'G6', 'G9', 'G12',
                                                                                 'H3', 'H6', 'H9', 'H12')  ~ 'FEC_Eim_Ct3'),
                                      Tm_sum = rowSums(eim_results[,22:25], na.rm = T),
                                      Tm_index = case_when(Well.Position %in% c('A1', 'A4', 'A7', 'A10', 
                                                                                'B1', 'B4', 'B7', 'B10', 
                                                                                'C1', 'C4', 'C7', 'C10', 
                                                                                'D1', 'D4', 'D7', 'D10',
                                                                                'E1', 'E4', 'E7', 'E10',
                                                                                'F1', 'F4', 'F7', 'F10',
                                                                                'G1', 'G4', 'G7', 'G10',
                                                                                'H1', 'H4', 'H7', 'H10') ~ 'Tm_sum1', 
                                                           Well.Position %in% c('A2', 'A5', 'A8', 'A11', 
                                                                                'B2', 'B5', 'B8', 'B11', 
                                                                                'C2', 'C5', 'C8', 'C11', 
                                                                                'D2', 'D5', 'D8', 'D11',
                                                                                'E2', 'E5', 'E8', 'E11',
                                                                                'F2', 'F5', 'F8', 'F11',
                                                                                'G2', 'G5', 'G8', 'G11',
                                                                                'H2', 'H5', 'H8', 'H11')  ~ 'Tm_sum2', 
                                                           Well.Position %in% c('A3', 'A6', 'A9', 'A12', 
                                                                                'B3', 'B6', 'B9', 'B12', 
                                                                                'C3', 'C6', 'C9', 'C12', 
                                                                                'D3', 'D6', 'D9', 'D12',
                                                                                'E3', 'E6', 'E9', 'E12',
                                                                                'F3', 'F6', 'F9', 'F12',
                                                                                'G3', 'G6', 'G9', 'G12',
                                                                                'H3', 'H6', 'H9', 'H12')  ~ 'Tm_sum3'),
                                      MC_eval = case_when(Tm_sum >= 79 ~ F, Tm_sum < 79 ~ T),
                                      MC_index = case_when(Well.Position %in% c('A1', 'A4', 'A7', 'A10', 
                                                                                'B1', 'B4', 'B7', 'B10', 
                                                                                'C1', 'C4', 'C7', 'C10', 
                                                                                'D1', 'D4', 'D7', 'D10',
                                                                                'E1', 'E4', 'E7', 'E10',
                                                                                'F1', 'F4', 'F7', 'F10',
                                                                                'G1', 'G4', 'G7', 'G10',
                                                                                'H1', 'H4', 'H7', 'H10') ~ 'MC.1', 
                                                           Well.Position %in% c('A2', 'A5', 'A8', 'A11', 
                                                                                'B2', 'B5', 'B8', 'B11', 
                                                                                'C2', 'C5', 'C8', 'C11', 
                                                                                'D2', 'D5', 'D8', 'D11',
                                                                                'E2', 'E5', 'E8', 'E11',
                                                                                'F2', 'F5', 'F8', 'F11',
                                                                                'G2', 'G5', 'G8', 'G11',
                                                                                'H2', 'H5', 'H8', 'H11')  ~ 'MC.2', 
                                                           Well.Position %in% c('A3', 'A6', 'A9', 'A12', 
                                                                                'B3', 'B6', 'B9', 'B12', 
                                                                                'C3', 'C6', 'C9', 'C12', 
                                                                                'D3', 'D6', 'D9', 'D12',
                                                                                'E3', 'E6', 'E9', 'E12',
                                                                                'F3', 'F6', 'F9', 'F12',
                                                                                'G3', 'G6', 'G9', 'G12',
                                                                                'H3', 'H6', 'H9', 'H12')  ~ 'MC.3'))


# pivot table to have all measurements per Mouse_ID in one row (not three)

## pivot by Ct_index
eim_results1 <- pivot_wider(eim_results, names_from = "Ct_index", values_from = "Cq") %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% 
  fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) %>% 
  select(-c(Tm_index, Tm_sum, MC_index, MC_eval))

# pivot by Tm_index
eim_results2 <- pivot_wider(eim_results, names_from = "Tm_index", values_from = "Tm_sum") %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% 
  fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) %>% 
  select(-c(Ct_index, MC_index, MC_eval, Tm1,Tm2,Tm3,Tm4, Cq))

# pivot by MC_index and introduce MC.Eimeria.FEC
eim_results3 <- pivot_wider(eim_results, names_from = "MC_index", values_from = "MC_eval") %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% 
  fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) %>% 
  select(-c(Ct_index,Tm_index,Tm_sum, Tm1,Tm2,Tm3,Tm4, Cq)) %>% 
  mutate(MC.Eimeria.FEC = case_when((MC.1 == T & MC.2 == T & MC.3 == T) ~ T,
                                    (MC.1 == F & MC.2 == T & MC.3 == T) | (MC.1 == T & MC.2 == F & MC.3 == T) | (MC.1 == T & MC.2 == T & MC.3 == F) ~ T,
                                    (MC.1 == F & MC.2 == F & MC.3 == T) | (MC.1 == T & MC.2 == F & MC.3 == F) | (MC.1 == F & MC.2 == T & MC.3 == F) ~ F,
                                    (MC.1 == F & MC.2 == F & MC.3 == F) ~ F))

# join the pivots
eim_results <- full_join(eim_results1, eim_results2)
eim_results <- full_join(eim_results, eim_results3)

# make Ct cols numeric 
# (There was the term "Undetermined" for samples without Ct measurements,
# that caused the column to be of type character)
eim_results$FEC_Eim_Ct1 <- as.numeric(eim_results$FEC_Eim_Ct1);
eim_results$FEC_Eim_Ct2 <- as.numeric(eim_results$FEC_Eim_Ct2);
eim_results$FEC_Eim_Ct3 <- as.numeric(eim_results$FEC_Eim_Ct3)

# (Un)Select Columns for final Data Product:
eim_results <- eim_results %>% 
  select(-c(Well, Well.Position, Omit, Quencher, Curve.Quality, Result.Quality.Issues, 
            Auto.Threshold)) %>%
  arrange(Mouse_ID)

# To check the efficiency of my MC.Eimeria.FEC column identification of true infections, 
# I compared a manually evaluated column (looking at the curve and giving T/F 
# based on the Melting Temp curve) to an automatically evaluated column as introduced above (MC.Eimeria.FEC)
      #MC.evaluation <- full_join(eim_results, weights_22) %>% select(Mouse_ID,Tm_sum1, Tm_sum2, Tm_sum3,MC.Eimeria.FEC, MC.Eimeria.FEC.vis)
      #MC.evaluation <- MC.evaluation %>% mutate(eval = MC.Eimeria.FEC == MC.Eimeria.FEC.vis)
      #write.csv(MC.evaluation, "~/Documents/Github/Crypto/HZ22/MC.evaluation.csv")


#write the data frame in a csv file --> you may need to use your local path
#write.csv(FEC_18_19_21_22, "./data_input/FEC_EqPCR_DataProduct_HZ18_HZ22.csv", row.names=FALSE)
write.csv(eim_results, "~/Documents/Github/Mouse_Eimeria_Field/data_input/FEC_EqPCR_input_data.csv", row.names=FALSE)


################################################################################
####  Crypto qPCR FEC  #########################################################
################################################################################

#creating script to read all the results file produced by the program QuantStudio 1
#the result files are located at: GitHub/Mouse_Eimeria_Field/data_input/qPCR_FEC_22/results

#create a list out of the names of each excel file (maybe you have to do that locally)
#crypto_plates <- list.files(path = "~/Documents/Github/Mouse_Eimeria_Field/data_input/qPCR_FEC_22/results/crypto_plates/")
#write.csv(crypto_plates, "~/Documents/Github/Mouse_Eimeria_Field/data_input/qPCR_FEC_22/results/crypto_plates.csv")
crypto_plates <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/qPCR_FEC_22/results/crypto_plates.csv")
crypto_list <- as.list(crypto_plates$x)

#write a function to specify how to read the qPCR files
read_cryp_qPCR_file <- function(x) {
  
  df1 <- xlsx::read.xlsx(x, sheetIndex = 1)
  #get the file name of the file
  #the name of the file is the second column of this file
  #to get that name we can start by selecting this column
  #and then getting the name out of it
  filename <- colnames(df1[2])
  #remove unecessary rows of the data frame
  #everything before actual data
  #different in Crypto (22) vs. Eimeria (23)!
  df1 <- df1 %>%
    filter(!row_number() %in% c(1:22))
  #change the column names to the names of the first row
  colnames(df1) <- df1[1, ]
  #Now remove the first row
  df1 <- df1 %>% filter(!row_number() %in% 1)
  #add a column with the name of the plate 
  df1 <- df1 %>% mutate(plate = filename)
  
}

#apply the function you created in the last step to each of the elements (names of files)
#of the list NT_list 
#in this way you are creating a list of data frames from the result files produced by the 
#program of the machine
setwd("/Users/finnlo/Documents/Github/Mouse_Eimeria_Field/data_input/qPCR_FEC_22/results/crypto_plates/")
results <- lapply(crypto_list, read_cryp_qPCR_file)

#show the data frame consisting of each result data file
cryp_results <- Reduce(rbind, results)

#remove duplicates (eliminates some NTCs)
cryp_results <- unique(cryp_results)

################################################################################
#### PIVOT #####################################################################
################################################################################

# correct colnames
colnames(cryp_results)[colnames(cryp_results)%in%"Sample"]  <- "Mouse_ID"
colnames(cryp_results)[colnames(cryp_results)%in%"Cq Mean"] <- "FEC_Crypto_Ct" # specific to this assay, to avoid confusion with other qPCRs

# correct colname for easier selecting
cryp_results <- dplyr::rename(cryp_results, Well_Position = "Well Position",
                              Curve_Quality = "Curve Quality",
                              Result_Quality_Issues = "Result Quality Issues",
                              Auto_Threshold = "Auto Threshold")

# correct type of some cols to numeric 
cryp_results[,c(1, 10, 13:16, 18, 20:21)] <- lapply(cryp_results[, c(1, 10, 13:16, 18, 20:21)], as.numeric)

# correct row input of Cq --> undetermined to 0 
# a 0 means that the sample was run, but nothing was measured, 
# whereas NA would mean the absence of the sample..
# BEWARE: For oocyst prediction this 0 will lead to a super high prediction, must be corrected to 0 again!
cryp_results <- cryp_results %>% mutate(FEC_Crypto_Ct = case_when(is.na(FEC_Crypto_Ct) ~ 0,
                                                                  !is.na(FEC_Crypto_Ct) ~ FEC_Crypto_Ct))
cryp_results$FEC_Crypto_Ct <- round(cryp_results$FEC_Crypto_Ct, 2)

# specify Ct measurements from well plate position
cryp_results <- cryp_results %>% mutate(Ct_index = case_when (Well_Position %in% c('A1', 'A4', 'A7', 'A10', 
                                                                                   'B1', 'B4', 'B7', 'B10', 
                                                                                   'C1', 'C4', 'C7', 'C10', 
                                                                                   'D1', 'D4', 'D7', 'D10',
                                                                                   'E1', 'E4', 'E7', 'E10',
                                                                                   'F1', 'F4', 'F7', 'F10',
                                                                                   'G1', 'G4', 'G7', 'G10',
                                                                                   'H1', 'H4', 'H7', 'H10') ~ 'FEC_Crypto_Ct1', 
                                                              Well_Position %in% c('A2', 'A5', 'A8', 'A11', 
                                                                                   'B2', 'B5', 'B8', 'B11', 
                                                                                   'C2', 'C5', 'C8', 'C11', 
                                                                                   'D2', 'D5', 'D8', 'D11',
                                                                                   'E2', 'E5', 'E8', 'E11',
                                                                                   'F2', 'F5', 'F8', 'F11',
                                                                                   'G2', 'G5', 'G8', 'G11',
                                                                                   'H2', 'H5', 'H8', 'H11')  ~ 'FEC_Crypto_Ct2', 
                                                              Well_Position %in% c('A3', 'A6', 'A9', 'A12', 
                                                                                   'B3', 'B6', 'B9', 'B12', 
                                                                                   'C3', 'C6', 'C9', 'C12', 
                                                                                   'D3', 'D6', 'D9', 'D12',
                                                                                   'E3', 'E6', 'E9', 'E12',
                                                                                   'F3', 'F6', 'F9', 'F12',
                                                                                   'G3', 'G6', 'G9', 'G12',
                                                                                   'H3', 'H6', 'H9', 'H12')  ~ 'FEC_Crypto_Ct3')) 


# pivot table to have all measurements per Mouse_ID in one row (not three)
cryp_results <- pivot_wider(cryp_results, names_from = "Ct_index", values_from = "Cq") %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% 
  fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T)

# make new cols numeric
cryp_results$FEC_Crypto_Ct1 <- as.numeric(cryp_results$FEC_Crypto_Ct1);
cryp_results$FEC_Crypto_Ct2 <- as.numeric(cryp_results$FEC_Crypto_Ct2);
cryp_results$FEC_Crypto_Ct3 <- as.numeric(cryp_results$FEC_Crypto_Ct3)

# (Un)Select Columns for final Data Product:
cryp_results <- cryp_results %>% 
  dplyr::select(-c(Well, Well_Position, Omit, Quencher, Curve_Quality, Result_Quality_Issues, Auto_Threshold))


cryp_results <- cryp_results %>% filter(Task == "UNKNOWN")
  
#write the data frame in a csv file 
#write.csv(cryp_results, "./data_input/FEC_CqPCR_DataProduct_HZ22.csv", row.names=FALSE)
write.csv(cryp_results, "~/Documents/Github/Mouse_Eimeria_Field/data_input/FEC_CqPCR_input_data.csv", row.names=FALSE)

