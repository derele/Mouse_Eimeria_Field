# Load libraries
library(tidyr)
library(dplyr)

#setwd("~/GitHub/Mouse_Eimeria_Field/")

# read the file that you need 
SOTA <- read.csv("data_products/SOTA_Data_Product.csv")


# ~20 extracted DNA samples of parents and hybrids 
# spanning the index  
# If no extracted DNA (from tissue, e.g. caecum or illeum) 
# liver or any other tissue would also be okay.
# sent them on the last week of November.
# 2019 mice would probably be easiest.

So_se <- SOTA %>%
  dplyr::filter(Year == "2019")

So_se$HI

# histogram to visualize samples
hist(So_se$HI)

# let's see what variables we have
glimpse(So_se)

#lets select the relevant columns to print out for sample hunting
So_se <- So_se %>% 
  dplyr::select(Mouse_ID, Year, HI)
