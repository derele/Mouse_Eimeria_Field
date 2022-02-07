library(xlsx)
library(tidyr)
library(dplyr)
library(janitor)

#read.cs
results <- read.xlsx("qPCR_eimeria_field_faeces_07022022_Results_20220207 123104.xlsx", sheetIndex = 1)
#remove unnecessary rows 
results <- results %>%
  filter(!row_number() %in% c(1:23))

colnames(results) <- results[1, ]

#Now remove the first row
results <- results %>% filter(!row_number() %in% 1)
