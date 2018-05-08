##Julia Analysis Flotation vs PCR

# Set folder where is the file of interest
setwd("../raw_data/Eimeria_detection/")

# read the file and give it a name 
pcr.data <- read.csv("Inventory_contents_all.csv")
 
#eliminate all the information that you won't need (optional)
pcr.data <- pcr.data[-c(5:11,14, 16, 18, 20,21)]

#change columns names (optional)
colnames(pcr.data) <- c("Year", "Transect", "Code", "Mouse_ID", "Flot", "Ap5", "n18S_Seq", 
                                                         "COI_Seq", "ORF470_Seq")
#create a new variable with date from one specific year
results.2017 <- subset(pcr.data, pcr.data$Year == "2017")
View(results.2017)
na.omit(results.2017$Ap5)

table(results.2017$Flot, useNA = "always")

finalData <- subset(x = results.2017,
                    subset = (results.2017$Flot %in% c("positive", "negative")) & !is.na(results.2017$Ap5))
                                         
length(which(finalData$Ap5 %in% "positive" &
               c(finalData$n18S_Seq == "positive" | finalData$COI_Seq == "positive" | finalData$ORF470_Seq == "positive")))
