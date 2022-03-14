library(xlsx)
library(tidyr)
library(dplyr)
library(janitor)
library(readr)

#creating script to read all the results file produced by the programm QuantStudio 1
#the result files are located at: GitHub/Mouse_Eimeria_Field/data_input/qPCR_CEWE_21/Results_Files/

#change back to the github repository (Change to where your github repository is located)
setwd("~/GitHub/")

#1 change the current working directory to the location where the qPCR data files are
setwd("Mouse_Eimeria_Field/data_input/qPCR_CEWE_21/Results_Files/")

#2. Create a list of names, out of the files in this repository
#list those files
list_CEWE <- as.list(list.files()) #now we have a list of the data frame names

#We need to create a list out of the names 
list_names <- as.vector(unlist(list_CEWE))

#change back to the github repository (Change to where your github repository is located)
setwd("~/GitHub/")

#change back to the github repository (Change to where your github repository is located)
setwd("Mouse_Eimeria_Field/data_input/qPCR_CEWE_21/Results_Files/")

#write a function to specify how to read the qPCR files
read_qPCR_file <- function(x) {
  
  df1 <- read.xlsx(x, sheetIndex = 1)
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
df_results <- unique(df_results)

#change your working directory
setwd("~/GitHub/Mouse_Eimeria_Field")

#write the data frame in a csv file 
write.csv(df_results, "data_products/qPCR_CEWE_eimeria_21_22", row.names=FALSE)
