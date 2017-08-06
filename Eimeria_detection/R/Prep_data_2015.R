# Import data
data2015 <- read.csv("../data_raw/Results_flotation_and_PCR_2015.csv")

# how many mice?
length(unique(data2015$Sample_ID))

#########################
## Export :
dataexport <- data.frame(Mouse_ID = data2015$Sample_ID,
                         OPG = data2015$OPG_.oocysts.ml.g_of_faeces.,
                         Year = 2015)

write.csv(x = dataexport, file = "../data_clean/Results_flotation_and_PCR_2015_CLEAN_temp.csv", row.names = FALSE)
