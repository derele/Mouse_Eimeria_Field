## Alice September 2017

## Import files to cross check
dissection2017 <- read.csv("../raw_data/HZ17_September_Mice_Dissection")
trapping2017 <- read.csv("../raw_data/HZ17_Mice_Trap.csv")

## Some kind of table adress/mouseid/data of sampling


# table : c(adress, code, long, lat, date)
 
 # print id_mouse
 
  # replace manually after check
 



## Extract individual coordinates from trapping
individual_mice_trapped <- unlist(strsplit(as.character(trapping2017$Mice_IDs), " "))
individual_mice_trapped <- gsub(pattern = " ", replacement = "", x = individual_mice_trapped)

individual_other_trapped <- unlist(strsplit(as.character(trapping2017$Other_samples), " "))
individual_other_trapped <- gsub(pattern = " ", replacement = "", x = individual_other_trapped)

link_table <- data.frame(Rodent_ID = c(individual_mice_trapped, individual_other_trapped),
                                        Longitude = NA, Latitude = NA)
for(i in 1:nrow(link_table)){
  if(trapping2017[grep(link_table$Rodent_ID[i], trapping2017$Mice_IDs),]$Longitude == TRUE){
    link_table$Longitude[i] <- trapping2017[grep(link_table$Rodent_ID[i], trapping2017$Mice_IDs),]$Longitude
  }
}


link_table$Longitude[i] 
link_table$Rodent_ID %in% trapping2017$Mice_IDs
  
  print(c(
  link_table$Latitude[i] ,
  trapping2017[grep(link_table$Rodent_ID[i], trapping2017$Mice_IDs),]$Latitude)
}


## Same for fishing in Ratfeces table...

trapping2017$Other_samples

## Cross check matching coordinates and Mouse_ID
dissection2017$Latitude %in% trapping2017$Latitude

dissection2017$Longitude %in% trapping2017$Longitude

