library(xlsx)
library(tidyr)
library(dplyr)
library(janitor)

#creating script to read all the results file produced by the programm QuantStudio 1
#the result files are located at: GitHub/Mouse_Eimeria_Field/data_input/qPCR_faeces/Results


#1 change location to where the qPCR data files are
setwd("/localstorage/fay/GitHub/Mouse_Eimeria_Field/data_input/qPCR_faeces/Results/Results_files")

#read the table whith the names of each file, call it NT (= Name table)
NT <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/qPCR_faeces/Results/Filenames_qpcr_results.csv")


NT_list <- as.list(NT$qPCR.Results.file.names)

#write a function to specify how to read the qPCR files
read_qPCR_file <- function(x) {
  setwd("~/GitHub/Mouse_Eimeria_Field/data_input/qPCR_faeces/Results/Results_files")
  df1 <- read.xlsx(x, sheetIndex = 1)
  #remove unecessary columns of the data frame
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

#change your working directory
setwd("~/GitHub/Mouse_Eimeria_Field")

#write the data frame in a csv file 
write.csv(df_results, "data_products/qPCR_faeces_2022", row.names=FALSE)
