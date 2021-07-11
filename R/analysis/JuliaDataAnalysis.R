##Julia Analysis Flotation vs PCR

# Set folder where is the file of interest
setwd("../raw_data/Eimeria_detection/")
setwd("/home/victor/Dokumente/Git_projects/Mouse_Eimeria_Databasing/raw_data/Eimeria_detection/")
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

write.csv(finalData, "/home/victor/Dokumente/Supervision/Julia/Eimeria_Detection_2017(final_data).csv")

#Total number of Ap5 positive                                         
length(which(finalData$Ap5 %in% "positive"))

#Total number of TRUE Ap5 positive 
length(which(finalData$Ap5 %in% "positive" &
               c(finalData$n18S_Seq == "positive" | finalData$COI_Seq == "positive" | finalData$ORF470_Seq == "positive")))

#Total number of Flotation positive
length(which(finalData$Flot %in% "positive"))

#Total number of TRUE Flotation positive 
length(which(finalData$Flot %in% "positive" & 
               c(finalData$n18S_Seq == "positive" | finalData$COI_Seq == "positive" | finalData$ORF470_Seq == "positive")))


##Venn diagram 
##Load libraries required 
library("VennDiagram")
library(grid)
library(gridExtra)

##Plot
grid.newpage()
draw.triple.venn(area1 = nrow(subset(finalData, n18S_Seq == "positive")), area2 = nrow(subset(finalData, COI_Seq == "positive")), 
                 area3 = nrow(subset(finalData, ORF470_Seq == "positive")), n12 = nrow(subset(finalData, n18S_Seq == "positive"&COI_Seq== "positive")), 
                 n23 = nrow(subset(finalData, COI_Seq== "positive"&ORF470_Seq =="positive")), n13 = nrow(subset(finalData, n18S_Seq == "positive"&ORF470_Seq == "positive")), 
                 n123 = nrow(subset(finalData, n18S_Seq == "positive"&ORF470_Seq == "positive"&COI_Seq=="positive")), category = c("18S", "COI", "ORF470"), 
                 lty = rep(1,3), col = c("dodgerblue4", "firebrick3", "darkgreen"), lwd = rep(2,3),
                 fill = c("dodgerblue4", "firebrick3", "darkgreen"), alpha = c(0.3, 0.3, 0.3), cex =2, cat.cex = 2.5, cat.default.pos = 'outer', 
                 cat.col = c("dodgerblue4", "firebrick3", "darkgreen"))

#Plot
grid.newpage()
draw.pairwise.venn(area1= length(which(finalData$Ap5 %in% "positive")), area2= length(which(finalData$Flot %in% "positive")), cross.area = length(which(finalData$Flot %in% "positive" & finalData$Ap5 %in% "positive")))
                   

#grid.newpage()
#draw.pairwise.venn(22, 20, 11, category = c("Dog People", "Cat People"), lty = rep("blank", 
#                                                                                  2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 
#                                                                                                                                                        0), cat.dist = rep(0.025, 2), scaled = FALSE)
#https://rstudio-pubs-static.s3.amazonaws.com/13301_6641d73cfac741a59c0a851feb99e98b.html