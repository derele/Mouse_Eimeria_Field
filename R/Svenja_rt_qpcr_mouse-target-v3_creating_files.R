library(ggplot2)
library(readxl)
library(dplyr)
#library(xlsx)
library(openxlsx)

folder<-"D:/UNI/DATA/RT_qPCR/definitivo/plate13" #selecting the plates one by one!
setwd(folder)
getwd()

plate13 <- read_excel("rt_qpcr_plate13.xls")


loaded_table<-plate13
df_loaded_table<- data.frame(loaded_table$Name, loaded_table$`Ct SYBR`, loaded_table$`Target SYBR`)
clean_table <- df_loaded_table[-c(85:101),]

#1. create the mean of all targets and housekeeping genes per ID
x<-(1:(nrow(clean_table)/3))#creating a list, where each ID appears just once(triplicats/3=1)

for (i in 1:(nrow(clean_table)/3)){
  ct1<- clean_table[i*3-2, ]$loaded_table..Ct.SYBR. #define at which position c1... is
  ct2<- clean_table[i*3-1, ]$loaded_table..Ct.SYBR.
  ct3<- clean_table[i*3, ]$loaded_table..Ct.SYBR.
  ct1<-as.numeric(levels(ct1)[ct1])
  if (is.na(ct1)){
    ct1 <- 1000
  }
  ct2<- as.numeric(levels(ct2)[ct2])
  if(is.na(ct2)){
    ct2 <- 1000
  }
  ct3<-as.numeric(levels(ct3)[ct3])
  if(is.na(ct3)){
    ct3 <- 1000
  }
  acc_range <- 1.5 #I accept no less than two ct's in acc_range
  
  if (abs(ct1-ct2)<= acc_range & abs(ct1-ct3) <= acc_range){
    
    x[i] <- (ct1+ct2+ct3)/3
  }else if (abs(ct1-ct2)<= acc_range & abs(ct1-ct3) > acc_range){
    
    x[i] <- (ct1+ct2)/2
  }else if (abs(ct1-ct2)> acc_range & abs(ct1-ct3) <= acc_range){
    
    x[i] <- (ct1+ct3)/2
  }else if (abs(ct2-ct3)<= acc_range & abs(ct1-ct2) > acc_range){
    
    x[i] <- (ct2+ct3)/2
  }else{
    
    x[i] <- 1000
  }
  
  
} # now, there are all accepted means stored in "x"

#2. create a dataframe with ID's, targets, means 
#ID's
liste_ID <- list()
y <- 1
for (i in clean_table$loaded_table.Name){ #i is a string now
  
  y <- y+1#y is now counting 
  if (y > 3) {
    print(i)  # every third name is printed
    liste_ID <- c(liste_ID,i) #every third name is stored in "liste"
    
    y <-1 #y ist 1 again
    
    
  }
}
mydata_ID <- do.call(rbind.data.frame, liste_ID) #is a df now

#Targets
liste_targets <- list()
y <- 1
for (i in clean_table$loaded_table..Target.SYBR.){ #i is a string now
  
  y <- y+1#y is now counting 
  if (y > 3) {
    print(i)  # every third target is printed
    liste_targets <- c(liste_targets,i) #every third target is stored in "liste"
    
    y <-1 #y ist 1 again
    
    
  }
}
mydata_targets <- do.call(rbind.data.frame, liste_targets) # is a df now
#create "means"
means <- dplyr::mutate(mydata_ID, x)
names(means) <- c("Mouse_ID", "deltaCt")
names(mydata_targets) #<- c("Mouse_ID", "deltaCt")
means <- dplyr::mutate(mydata_ID,liste_targets, x)
names(means) <- c("Name","Target SYBR", "Ct mean")
#create subsets
GAPDH <- subset(means, means$`Target SYBR` == "GAPDH")
beta_actin <- subset(means, means$`Target SYBR` == "beta-Actin")
targets <- subset(means, means$`Target SYBR` != "GAPDH")
targets <- subset(targets, targets$`Target SYBR` != "beta-Actin")

#3. calculate delta ct = GAPDH-target

xy<-list()
x <- 1
for (i in 1:20){
  print(x)
  ct1 <-GAPDH[ x, ]$`Ct mean`
  ct2<-targets[i, ]$`Ct mean`
  ct1<-as.numeric(ct1)
  ct2<-as.numeric(ct2)
  x <- x+1
  xy<- c(xy,ct1 -ct2)
  if (x >4) {x <-1}
}

targets <- dplyr::mutate(targets, xy)
names(targets) <- c("Name","Target SYBR", "Ct mean","delta Ct GAPDH")
#4. calculate delta ct = beta actin-target

xy<-list()
x <- 1
for (i in 1:20){
  print(x)
  ct1 <-beta_actin[ x, ]$`Ct mean`
  ct2<-targets[i, ]$`Ct mean`
  ct1<-as.numeric(ct1)
  ct2<-as.numeric(ct2)
  x <- x+1
  xy<- c(xy,ct1 -ct2)
  if (x >4) {x <-1}
}

targets <- dplyr::mutate(targets, xy)
names(targets) <- c("Name","Target SYBR", "Ct mean","delta Ct GAPDH", "delta Ct beta actin")

#creating an index: (GAPDH+betaActin)/2
index <- list()
b <- 1
for (i in 1:20){
 print(i) 
  b <- (targets$`delta Ct GAPDH`[[i]]+targets$`delta Ct beta actin`[[i]])/2
 # b <- b+1
 # st <- c(st, b)
  print(b)
  index <- c(index, b)
  
}

targets <- dplyr::mutate(targets, index)


#5. export the file. Further proessing in rtqpcr.R
tmp = paste(folder, "rtResults.xlsx", sep = "/")  
write.xlsx(targets, tmp)









