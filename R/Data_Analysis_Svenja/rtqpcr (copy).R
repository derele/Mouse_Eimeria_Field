library(ggplot2)
library(readxl)
library(dplyr)
library(openxlsx)

folder<-"D:/UNI/DATA/RT_qPCR/definitivo"
setwd(folder)
getwd()
source("../../../R codes/continious_axis_delta_ct_with_BLAST.R")
#files created in: rt_qpcr_mouse-target-v3_creating_files.R
plate1 <- read_excel("plate1/rtResults.xlsx")
plate1 <- plate1[!grepl("CEWE_AA_531", plate1$Name),]
plate2<- read_excel("plate2/rtResults.xlsx")
plate3 <- read_excel("plate3/rtResults.xlsx")
plate3 <- plate3[!grepl("CEWE_AA_534", plate3$Name),]
plate4 <- read_excel("plate4/rtResults.xlsx")
plate5 <- read_excel("plate5/rtResults.xlsx")
plate6 <- read_excel("plate6/rtResults.xlsx")
plate7 <- read_excel("plate7/rtResults.xlsx")
plate8 <- read_excel("plate8/rtResults.xlsx")
plate8 <- plate8[!grepl("CEWE_AA_569", plate8$Name),]
plate8 <- plate8[!grepl("CEWE_AA_635", plate8$Name),]
plate8 <- plate8[!grepl("CEWE_AA_561", plate8$Name),]
plate9 <- read_excel("plate9/rtResults.xlsx")
plate10 <- read_excel("plate10/rtResults.xlsx")
plate10 <- plate10[!grepl("CEWE_AA_591", plate10$Name),]
plate11 <- read_excel("plate11/rtResults.xlsx")
plate11 <- plate11[!grepl("CEWE_AA_636", plate11$Name),]
plate12 <- read_excel("plate12/rtResults.xlsx")
plate12 <- plate12[!grepl("CEWE_AA_662", plate12$Name),]
plate13 <- read_excel("plate13/rtResults.xlsx")
plate14 <- read_excel("plate14/rtResults.xlsx")


total <- rbind(plate1,plate2,plate3,plate4,plate5,plate6,plate7,
               plate8,plate9,plate10,plate11,plate12,plate13,plate14)
total$`delta Ct GAPDH` <-as.numeric(total$`delta Ct GAPDH`)
total$`delta Ct beta actin` <- as.numeric(total$`delta Ct beta actin`)
total$index <- as.numeric(total$index)
result <- read_excel("../../../DATA/Results.xlsx")
result$Mouse_ID <- gsub(pattern = "AA_0", replacement = "CEWE_AA_", x = result$Mouse_ID)




total$eimeria_cecum     <- "FALSE" #add a coloumn with "eimeria_cecum" and store "FALSE" in each row

for (i in 1:dim(total)[1]){  #for loop to add the right eimeria_cecum from "result" in "total"
  
  for (j in 1:dim(result)[1]){
    
    if (result[j,1]==total[i,1]){
      total[i,7] <- result[j,4]
      return
    }
  }
}
total$status <- "FALSE"
for (i in 1:dim(total)[1]){  #for loop to add the right status from "result" in "total"
  
  for (j in 1:dim(result)[1]){
    
    if (result[j,1]==total[i,1]){
      total[i,8] <- result[j,9]
      return
    }
  }
}
total$worms <- "FALSE"
for (i in 1:dim(total)[1]){  #for loop to add the right worm status from "result" in "total"
  
  for (j in 1:dim(result)[1]){
    
    if (result[j,1]==total[i,1]){
      total[i,9] <- result[j,5]
      return
    }
  }
}
total$crypto <- "FALSE"
for (i in 1:dim(total)[1]){  #for loop to add the right crypto status from "result" in "total"
  
  for (j in 1:dim(result)[1]){
    
    if (result[j,1]==total[i,1]){
      total[i,10] <- result[j,3]
      return
    }
  }
}
total$ct_eimeria <- "FALSE"
delta_ct_eimeria$Name <- as.character(delta_ct_eimeria$Name)
for (i in 1:dim(total)[1]){  #for loop to add the ct from eimeria from "result" in "total"
  
  for (j in 1:dim(delta_ct_eimeria)[1]){
    
    if (delta_ct_eimeria[j,1]==total[i,1]){
      total[i,11] <- delta_ct_eimeria[j,2]
      return
    }
  }
}
total$eimeria_subspecies <- "FALSE"
for (i in 1:dim(total)[1]){  #for loop to add the eimeria subspecies from "result" in "total"
  
  for (j in 1:dim(result)[1]){
    
    if (result[j,1]==total[i,1]){
      total[i,12] <- result[j,10]
      return
    }
  }
}
total$ct_eimeria<- as.numeric(total$ct_eimeria)
####export total###########
tmp = paste(folder, "total.xlsx", sep = "/")  
write.xlsx(total, tmp)