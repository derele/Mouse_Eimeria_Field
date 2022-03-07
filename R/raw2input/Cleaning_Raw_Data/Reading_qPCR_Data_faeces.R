library(xlsx)
library(tidyr)
library(dplyr)
library(janitor)
library(readr)

#creating script to read all the results file produced by the programm QuantStudio 1
#the result files are located at: GitHub/Mouse_Eimeria_Field/data_input/qPCR_faeces/Results

#change back to the github repository (Change to where your github repository is located)
setwd("~/GitHub/")

#1 change the current working directory to the location where the qPCR data files are
setwd("Mouse_Eimeria_Field/data_input/qPCR_faeces/Results/Results_files/")

#read the table whith the names of each file, call it NT (= Name table)
NT <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/qPCR_faeces/Results/Filenames_qpcr_results.csv")

#create a list out of the names of each excel file
NT_list <- as.list(NT$qPCR.Results.file.names)

#change back to the github repository (Change to where your github repository is located)
setwd("~/GitHub/")

#change back to the github repository (Change to where your github repository is located)
setwd("Mouse_Eimeria_Field/data_input/qPCR_faeces/Results/Results_files")

#write a function to specify how to read the qPCR files
read_qPCR_file <- function(x) {
  
  df1 <- read.xlsx(x, sheetIndex = 1)
  #remove unecessary rows of the data frame
  #everything before actual data
  df1 <- df1 %>%
    filter(!row_number() %in% c(1:23))
  #change the column names to the names of the first row
  colnames(df1) <- df1[1, ]
  #Now remove the first row
  df1 <- df1 %>% filter(!row_number() %in% 1)
  
}


#apply the function you created in the last step to each of the elements (names of files)
#of the list NT_list 
#in this way you are creating a list of data frames from the result files produced by the 
#program of the machine
list_results <- lapply(NT_list, read_qPCR_file)

#show the data frame consisting of each result data file
df_results <- Reduce(rbind, list_results)

df_results <- unique(df_results)

#change your working directory
setwd("~/GitHub/Mouse_Eimeria_Field")

#write the data frame in a csv file 
write.csv(df_results, "data_products/qPCR_faeces_field_eimeria_22", row.names=FALSE)

