library(readxl)
library(dplyr)
library(tidyr)
library(stringr)

HZ21_oocysts <- read.csv("Mouse_Eimeria_Field/data_input/HZ21_Oocysts.csv")

sota <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/SOTA_Data_Product.csv")

basics          <- c("Mouse_ID", "Address", "Sex", "Longitude", 
                     "Latitude", "Year", "HI", "HI_NLoci")

oocyst.cols     <- c("counter", "Feces_Weight", "c", "N_oocysts_sq1",
                     "N_oocysts_sq2", "N_oocysts_sq3",  "N_oocysts_sq4",
                     "N_oocysts_sq5", "N_oocysts_sq6", "N_oocysts_sq7",
                     "N_oocysts_sq8", "mean_neubauer", "PBS_dil_in_mL", 
                     "OPG", "Ncells")

#concatenate the two vectors
basics_ooc <- c(basics, oocyst.cols)

#select from sota the relevant columns
sota_selection <- sota %>% select(all_of(basics_ooc))
glimpse(sota_selection)

### Data Cleaning
## View the data
#check the type of each column
glimpse(HZ21_oocysts)

#mouse ids
HZ21_oocysts$Mouse_ID <- gsub("AA_", "AA_0", HZ21_oocysts$Mouse_ID)

# feces_g: 1. Column name should be Feces_Weight 
#2. type should be numeric 
HZ21_oocysts <- HZ21_oocysts %>% dplyr::rename(Feces_Weight = Feces_g)
HZ21_oocysts$Feces_Weight <- as.numeric(HZ21_oocysts$Feces_Weight)

#year should be Year
HZ21_oocysts <- HZ21_oocysts %>% dplyr::rename(Year = year)

#what about the dates?
#in sota dates are month year 
HZ21_oocysts$date_count[grep("11/*.", HZ21_oocysts$date_count)] <- "November2021"
HZ21_oocysts$date_count[grep("12/*.", HZ21_oocysts$date_count)] <- "December2021"

#date_count should be Date_count
HZ21_oocysts <- HZ21_oocysts %>% dplyr::rename(Date_count = date_count)

#column ncells is incorrectly labelled
typeof(HZ21_oocysts$NCells)
typeof(sota_selection$Ncells) #Ncells should be an integer

HZ21_oocysts <- HZ21_oocysts %>% dplyr::rename(Ncells = NCells)
HZ21_oocysts$Ncells <- as.integer((HZ21_oocysts$Ncells))

#what about the dilutions?
typeof(sota_selection$PBS_dil_in_mL) #type is double
typeof(HZ21_oocysts$dilution_ml) #type is integer

#rename to the correct name
HZ21_oocysts <- HZ21_oocysts %>% dplyr::rename(PBS_dil_in_mL = dilution_ml)

#change the type
HZ21_oocysts$PBS_dil_in_mL <- as.double((HZ21_oocysts$PBS_dil_in_mL))

#the opg sometimes contains "Value"
typeof(HZ21_oocysts$OPG) #is a character

#create new column for the oocysts per gram to replace the old
HZ21_oocysts <- HZ21_oocysts %>% 
  mutate(OPG = (((N_oocysts_sq1 + N_oocysts_sq2 + N_oocysts_sq3 + N_oocysts_sq4 
                    + N_oocysts_sq5 + N_oocysts_sq6 + N_oocysts_sq7 + N_oocysts_sq8) / Ncells) * 
                    (PBS_dil_in_mL / (0.0001 * Ncells))) 
         / Feces_Weight)


#create new column to calculate the mean neubauer
HZ21_oocysts <- HZ21_oocysts %>% 
  mutate(mean_neubauer = (N_oocysts_sq1 + N_oocysts_sq2 + N_oocysts_sq3 + N_oocysts_sq4 
                  + N_oocysts_sq5 + N_oocysts_sq6 + N_oocysts_sq7 + N_oocysts_sq8) / 8)

#remove unecessary column N_oocysts_all_squares
HZ21_oocysts <- HZ21_oocysts %>% select(-N_oocysts_all_sq)

#write the cleaned table in the repository

write.csv(HZ21_oocysts, "Mouse_Eimeria_Field/data_input/Eimeria_detection/HZ21_Oocysts_cleaned.csv", row.names=FALSE)

