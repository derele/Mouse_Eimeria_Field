library(tidyverse)
library(data.table)
library(visdat)

#setwd("/Users/FinnLo/Documents/Programming/R/HZ_SC_and_Raw_Data/Mouse_Eimeria_Field/data_input/")
#### Creating a STATE OF THE ART Data Product (SOTA Data Product) ##############
################################################################################


#### Select Columns ############################################################
basics          <- c("Mouse_ID", "Address", "Sex", "Longitude", 
                     "Latitude", "Year", "HI", "HI_NLoci")

gen.loci        <- c("mtBamH", "YNPAR", "X332", "X347", "X65", "Tsx", "Btk", "Syap1",
                     "Es1", "Gpd1", "Idh1", "Mpi", "Np", "Sod1", "Es1C", "Gpd1C",
                     "Idh1C", "MpiC", "NpC", "Sod1C", "HI_NLoci",
                     "HI", "Zfy2", "Y", "Region")

dissection.cols <- c("Body_Weight", "Body_Length", "Tail_Length", "Status", "Spleen", 
                     "Left_Testis", "Right_Testis", "Seminal_Vesicles_Weight", "Liver",
                     "Sperm", "Left_Epididymis", "Right_Epididymis", 
                     "Right_Ovarium_Weight", "Left_Ovarium_Weight",
                     "Left_Embryo", "Right_Embryo", "Fleas", "Ticks", "Ectoparasites_Logical", 
                     "Arrival", "Dissection_Date", "Trap_Date", "Host")

tissue.cols     <- c("SPL1", "SPL2", "ELFO", "LIV", "KI", "LUN", "SG",
                     "MES", "COWE", "COCE", "COCE2", "CEWE", "CECE", 
                     "ILWE", "SICE", "WEOH", "WFOR", "FEC")

initial.worms.cols   <- c("Aspiculuris","Syphacia_obvelata","Trichuris_muris", "Taenia_taeniformis", "Flea", "Ticks", "Fleas",
                             "Ectoparasites",    "Worms_presence", "Syphacia",
                             "Hymenolepis_diminiuta", "Hymenolepis_diminuta", "Taenia_martis",    "Heligmosomoides_polygurus" ,
                             "Heterakis_spumosa","Mastophorus_muris","Hymenolepis_microstoma", "Catenotaenia_pusilla",
                             "Cysticercus", "Hymenolepis", "Taenia", "Aspiculuris_Syphacia",  
                             "Trichuris", "Heterakis", "Mastophorus", "Ectoparasites_Logical", 
                             "Aspiculuris", "Catenotaenia_pusilla")

final.worms.cols <- c("Aspiculuris_sp", "Syphacia_sp", "Trichuris_muris", "Taenia_sp",
                      "Heterakis_sp", "Mastophorus_muris", "Hymenolepis_sp", "Catenotaenia_pusilla",
                      "Heligmosomoides_polygurus", "Worms_presence")

oocyst.cols     <- c("counter", "Feces_Weight", "Date_count", "N_oocysts_sq1",
                     "N_oocysts_sq2", "N_oocysts_sq3",  "N_oocysts_sq4",
                     "N_oocysts_sq5", "N_oocysts_sq6", "N_oocysts_sq7",
                     "N_oocysts_sq8", "mean_neubauer", "PBS_dil_in_mL", 
                     "OPG", "Ncells")

EqPCR.cols      <- c("delta_ct_ilwe_MminusE", "delta_ct_cewe_MminusE", "MC.Eimeria", "Ct.Eimeria", "Ct.Mus")

EimGeno.cols    <- c("n18S_Seq", "COI_Seq", "ORF470_Seq", "eimeriaSpecies")

Gene.Exp.cols   <- c("IFNy",        "IL.12",       "IRG6",        "CXCR3",           
                     "IL.6"    ,    "GBP2")

CellCount.cols <- c( "Treg",   "CD4", "Treg17",      "Th1",         "Th17",   
                     "CD8",    "Act_CD8",     "IFNy_CD4",    "IL17A_CD4",  
                     "IFNy_CD8")

Crypto_qPCR.cols <- c("Ct_mean", "Oocyst_Predict")


#### Work Flow #################################################################
    ## Since Alice already put together a pretty extensive data product for the years
    ## until 2017, we will use this as a starting point.
    ## Alice was very un-invasive, so she left the columns unchanged, even if they
    ## contain the same data point, but column names differed slightly (e.g. spelling)
    ## We will clean that up, narrow down the columns to a manageable size.
    
    ## Once there are no double columns left, we can continue to add new data 
    ## (that also complies with our standardized column names)
    ## and add new assays, ...

#### 1. Load Data #################################################################
Alice <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/MiceTable_fullEimeriaInfos_2014to2017.csv")
Alice$HI_NLoci <- gsub(pattern = "HI ", replacement = "", x = Alice$HI_NLoci)
Alice$HI_NLoci <- as.integer(Alice$HI_NLoci)
Alice$Mouse_ID <- gsub(pattern = "Sk3173", replacement = "SK_3173", x = Alice$Mouse_ID)
  ## remove "useless" mice
wsh <- c(paste0("AA_000", 1:9), paste0("AA_00", 10:46))
apd <- c("A_0001", "A_0002", "A_0003")
useless <- c(wsh, apd)
Alice <- Alice[!(Alice$Mouse_ID %in% useless),]

Jarda <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ19_GenotypingJarda.csv", na.strings=c(""," ","NA"))
setnames(Jarda, old = c("PIN", "X_Longit", "Y_Latit"), new = c("Mouse_ID", "Longitude", "Latitude"), skip_absent = T)
Jarda$Mouse_ID <- gsub(pattern = "SK", replacement = "SK_", x = Jarda$Mouse_ID)
Jarda$Mouse_ID <- gsub(pattern = "Sk3173", replacement = "SK_3173", x = Jarda$Mouse_ID)

  ## merge and remove duplicates
SOTA <- full_join(Alice, Jarda[colnames(Jarda) %in% c(basics, gen.loci)]) %>% select(!which(!rowSums(!is.na(Alice)))) %>% select(!which(!colSums(!is.na(Alice)))) %>%
  arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 


################################################################################
#### 2. Time for Column Correction #############################################

## Address
    ## In this case, there are two columns called Address and Locality that essentially
    ## contain exactly the same info, but they were named differently in different years
    ## we will combine this data in one column by using pivot_longer in a temporary data frame
    ## this allows us to keep all data per Mouse_ID, and afterwards we can join the data
    ## with the original "SOTA" data frame
    ## in this process we will also look for partial duplicates (sometimes there are 2 rows of the same Mouse_ID
    ## that couldn't be joined because some row content differs), we will compare the rows, fill in missing (NA)
    ## data and then remove the (now "full") duplicates.
    ## Afterwards, we deselect unneeded columns
MT_Address <- SOTA %>% select(Mouse_ID, Address, Locality) %>% 
  pivot_longer(names_to = "Temp", values_to = "Address", cols = c(Address, Locality)) %>% 
  arrange(Mouse_ID) %>% 
  group_by(Mouse_ID) %>% 
  fill(c(everything()), .direction = "downup") %>% 
  ungroup() %>% 
  select(Mouse_ID, Address) %>% 
  distinct(Mouse_ID, .keep_all = T) 
SOTA <- full_join(SOTA, MT_Address) %>% 
  select(-c(Locality)) %>% 
  arrange(Mouse_ID) %>% 
  group_by(Mouse_ID) %>% 
  fill(c(everything()), .direction = "downup") %>% 
  ungroup() %>% 
  distinct(Mouse_ID, .keep_all = T)
rm(MT_Address)

## Aspiculuris
colnames(SOTA)[colnames(SOTA)%in%"Aspiculuris_tetraptera"] <- "Aspiculuris"


## Body_Weight
MT_Body_Weight <- SOTA %>% select(Mouse_ID, BW, Body_weight)
MT_Body_Weight <- MT_Body_Weight %>% pivot_longer(names_to = "Temp", values_to = "Body_Weight", cols = c(Body_weight, BW)) %>% 
  arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% select(Mouse_ID, Body_Weight) %>% distinct(Mouse_ID, .keep_all = T) 
SOTA <- full_join(SOTA, MT_Body_Weight) %>% select(-c(Body_weight, BW))
rm(MT_Body_Weight)
SOTA <- SOTA%>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 


## Body_Length == "Body_length", "Body_length1" -----> correct one super low value 8.9 or something.. --> 8.9
MT_Body_Length <- SOTA %>% select(Mouse_ID, L, Body_length)
MT_Body_Length <- MT_Body_Length %>% pivot_longer(names_to = "Temp", values_to = "Body_Length", cols = c(L, Body_length)) %>% 
  arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% select(Mouse_ID, Body_Length) %>% distinct(Mouse_ID, .keep_all = T) 
SOTA <- full_join(SOTA, MT_Body_Length) %>% select(-c(Body_length))
rm(MT_Body_Length)
SOTA <- SOTA%>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 

## Dissection_Date
colnames(SOTA)[colnames(SOTA)%in%"Dissection"] <- "Dissection_Date"

## Embryos
colnames(SOTA)[colnames(SOTA)%in%"Embryo_left"] <- "Left_Embryo"
colnames(SOTA)[colnames(SOTA)%in%"Embryo_right"] <- "Right_Embryo"


## Ectoparasites == "Ectoparasites", "Ectoparasites_Logical"
SOTA$Ectoparasites_Logical <- as.logical(SOTA$Ectoparasites)
MT_Ectoparasites_Logical <- SOTA %>% select(Mouse_ID, Ectoparasites) %>% pivot_longer(names_to = "Temp", values_to = "Ectoparasites_Logical", cols = c(Ectoparasites)) %>% 
  arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% select(Mouse_ID, Ectoparasites_Logical) %>% distinct(Mouse_ID, .keep_all = T) 
MT_Ectoparasites_Logical$Ectoparasites_Logical <- as.logical(MT_Ectoparasites_Logical$Ectoparasites_Logical)
SOTA <- full_join(SOTA, MT_Ectoparasites_Logical) %>% select(-c(Ectoparasites))
rm(MT_Ectoparasites_Logical)
SOTA <- SOTA %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 


## Epididymis
## Left_Epididymis == "Left_epididymis", "Left.epididymis.weight"
SOTA$Left.epididymis.weight <- as.double(SOTA$Left.epididymis.weight)
MT_Left_Epididymis <- SOTA %>% select(Mouse_ID, Left.epididymis.weight) %>% pivot_longer(names_to = "Temp", values_to = "Left_Epididymis", cols = c(Left.epididymis.weight)) %>% 
  arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% select(Mouse_ID, Left_Epididymis) %>% distinct(Mouse_ID, .keep_all = T) 
SOTA <- full_join(SOTA, MT_Left_Epididymis) %>% select(-c(Left.epididymis.weight))
rm(MT_Left_Epididymis)
SOTA <- SOTA%>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 


## Feces_Weight == "Feces_weight", "Feces_g", Feces
SOTA$Feces_weight <- as.double(SOTA$Feces_weight)
MT_Feces_Weight <- SOTA %>% select(Mouse_ID, Feces_weight, Feces_g) %>% pivot_longer(names_to = "Temp",  values_to = "Feces_Weight", cols = c(Feces_weight, Feces_g)) %>%
  arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% select(Mouse_ID, Feces_Weight) %>% distinct(Mouse_ID, .keep_all = T) 
SOTA <- full_join(SOTA, MT_Feces_Weight) %>% select(-c(Feces_g, Feces_weight))
rm(MT_Feces_Weight)
SOTA <- SOTA%>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 


## Fleas == "Flea", "Fleas"
SOTA <- SOTA %>% mutate(Fleas_Count = ifelse(Flea %in% c("0", "1", "2", "3", "4", "5", "6", "9", "11", "12"), as.numeric(Flea),
                                             ifelse(NA)),
                        Fleas = case_when(Flea == "fleas" ~ T,
                                          Flea == "TRUE" ~ T,
                                          Flea == "TRUE (collected)" ~ T,
                                          Flea == "FALSE"~ F,
                                          Fleas_Count > 0 ~ T,
                                          Fleas_Count == 0 ~ F,
                                          Ectoparasites_Logical == T ~ T,
                                          Ectoparasites_Logical == F ~ F,
                                          Ectoparasites_Logical == T & is.na(Flea) ~ T)) %>%
                        
  select(-c(Flea, Fleas_Count))
SOTA <- SOTA%>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 

## Head_Taken == "Head_taken", "Head.taken."
MT_Head_Taken <- SOTA %>% select(Mouse_ID, Head.taken.) %>%  pivot_longer(names_to = "Temp", values_to = "Head_Taken", cols = c(Head.taken.)) %>% 
  arrange(Mouse_ID) %>% group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Head_Taken)
MT_Head_Taken <- MT_Head_Taken %>% distinct(Mouse_ID, .keep_all = T) 
SOTA <- full_join(SOTA, MT_Head_Taken) %>% select(-c(Head.taken.))
rm(MT_Head_Taken)


## Latitude == "Latitude", "longitude" ***** MIX UP WITH LONG *****
MT_Latitude <- SOTA %>% select(Mouse_ID, Latitude, longitude) %>% pivot_longer(names_to = "Temp",  values_to = "Latitude", cols = c(Latitude, longitude)) %>%
  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%   fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Latitude)
MT_Latitude <- MT_Latitude %>% distinct(Mouse_ID, .keep_all = T) 
SOTA <- full_join(SOTA, MT_Latitude) %>% select(-c(longitude))
rm(MT_Latitude)


## Longitude == "Longitude", "latitude" ***** MIX UP WITH LAT *****
MT_Longitude <- SOTA %>% select(Mouse_ID, Longitude, latitude)  %>% pivot_longer(names_to = "Temp", values_to = "Longitude", cols = c(Longitude, latitude)) %>% 
  arrange(Mouse_ID) %>% group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Longitude)
MT_Longitude <- MT_Longitude %>% distinct(Mouse_ID, .keep_all = T) 
## join
SOTA <- full_join(SOTA, MT_Longitude) %>% select(-c(latitude))
rm(MT_Longitude)


## Liver == "Liver", "Liver_mass"
MT_Liver <- SOTA %>% select(Mouse_ID, Liver_mass) %>% pivot_longer(names_to = "Temp",  values_to = "Liver", cols = c(Liver_mass)) %>%  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Liver)
MT_Liver <- MT_Liver %>% distinct(Mouse_ID, .keep_all = T) 
SOTA <- full_join(SOTA, MT_Liver) %>% select(-c(Liver_mass))
rm(MT_Liver)

## Notes == "comments", "Note", "Notes", ("Embryo")
MT_Notes <- SOTA %>% select(Mouse_ID, Note, Notes, comments)  %>% pivot_longer(names_to = "Temp",  values_to = "Notes", cols = c(Note, Notes, comments)) %>%  
  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Notes)
MT_Notes <- MT_Notes %>% distinct(Mouse_ID, .keep_all = T) 
SOTA <- full_join(SOTA, MT_Notes) %>% select(-c(Note, comments))
rm(MT_Notes)
SOTA <- SOTA %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 


## Ovaria ==
## Right_Ovarium_Weight == "Right_ovarium", Right.ovarium.weight"
MT_Right_Ovarium_Weight <- SOTA %>% select(Mouse_ID, Right.ovarium.weight) %>% pivot_longer(names_to = "Temp",  values_to = "Right_Ovarium_Weight", cols = c(Right.ovarium.weight)) %>%  
  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Right_Ovarium_Weight) %>% distinct(Mouse_ID, .keep_all = T) 
SOTA <- full_join(SOTA, MT_Right_Ovarium_Weight) %>% select(-c(Right.ovarium.weight))
rm(MT_Right_Ovarium_Weight)


## Left_Ovarium_Weight == "Left_ovarium", "Left.ovarium.weight"
MT_Left_Ovarium_Weight <- SOTA %>% select(Mouse_ID, Left.ovarium.weight) %>% pivot_longer(names_to = "Temp",  values_to = "Left_Ovarium_Weight", cols = c(Left.ovarium.weight)) %>%  
  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Left_Ovarium_Weight) %>% distinct(Mouse_ID, .keep_all = T) 
SOTA <- full_join(SOTA, MT_Left_Ovarium_Weight) %>% select(-c(Left.ovarium.weight))
rm(MT_Left_Ovarium_Weight)


## Region == "Region", "REGion"
MT_Region <- SOTA %>% select(Mouse_ID, Region, REGion) %>% pivot_longer(names_to = "Temp",  values_to = "Region", cols = c(Region, REGion)) %>%
  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Region) %>% distinct(Mouse_ID, .keep_all = T) 
SOTA <- full_join(SOTA, MT_Region) %>% select(-c(REGion))
rm(MT_Region)
SOTA <- SOTA %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 


## Seminal_Vesicles_Weight == SemVes, Seminal.vesicle.weight, Seminal_Vesicles_Weight
MT_Seminal_Vesicles_Weight <- SOTA %>% select(Mouse_ID, SemVes, Seminal.vesicle.weight) %>% pivot_longer(names_to = "Temp",  values_to = "Seminal_Vesicles_Weight", cols = c(SemVes, Seminal.vesicle.weight)) %>%
  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Seminal_Vesicles_Weight) %>% distinct(Mouse_ID, .keep_all = T) 
SOTA <- full_join(SOTA, MT_Seminal_Vesicles_Weight) %>% select(-c(SemVes, Seminal.vesicle.weight))
rm(MT_Seminal_Vesicles_Weight)


## Spleen == "Spleen", "Spleen_mass"
SOTA$Spleen_mass <- as.double(SOTA$Spleen_mass)
MT_Spleen <- SOTA %>% select(Mouse_ID, Spleen, Spleen_mass)  %>% pivot_longer(names_to = "Temp",  values_to = "Spleen",cols = c(Spleen, Spleen_mass)) %>% 
  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Spleen) %>% distinct(Mouse_ID, .keep_all = T) 
SOTA <- full_join(SOTA, MT_Spleen) %>% select(-c(Spleen_mass))
rm(MT_Spleen)
SOTA <- SOTA %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 


## Syphacia
colnames(SOTA)[colnames(SOTA)%in%"Syphacia_obvelata"] <- "Syphacia"



## Testis_mass Separation
## Wrong data input for "Testis_mass": instead of individual Left_Testis or 
## Right_Testis data, a combination of "Left_Testis/Right_Testis" was supplied
## Separation of that Column into 2 Columns, original "Testis"  col. discarded
SOTA_Sep <- SOTA %>% select(Mouse_ID, Testis_mass) %>% filter(Testis_mass != "NA")
SOTA_Sep <- setDT(SOTA_Sep)[, paste0("Testis_mass", 1:2) := tstrsplit(Testis_mass, "/")]
setnames(SOTA_Sep, old = c("Testis_mass1", "Testis_mass2"), new = c("Left_Testis1", "Right_Testis1"), skip_absent = T)
SOTA_Sep$Left_Testis1 <- as.double(SOTA_Sep$Left_Testis1)
SOTA_Sep$Right_Testis1 <- as.double(SOTA_Sep$Right_Testis1)
SOTA <- full_join(SOTA, SOTA_Sep) %>% select(-Testis_mass)
rm(SOTA_Sep)


    ## Left_Testis  == "Left_Testis1", "Left_Testis_mass"
MT_Left_Testis <- SOTA %>% select(Mouse_ID, Left_Testis1, Left_Testis_mass) %>% pivot_longer(names_to = "Temp",  values_to = "Left_Testis", cols = c(Left_Testis1, Left_Testis_mass)) %>% 
  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Left_Testis) %>% distinct(Mouse_ID, .keep_all = T) 
SOTA <- full_join(SOTA, MT_Left_Testis) %>% select(-c(Left_Testis_mass, Left_Testis1))
rm(MT_Left_Testis)


    ## Right_Testis == Right_Testis1", "Right_Testis_mass"
MT_Right_Testis <- SOTA %>% select(Mouse_ID, Right_Testis1, Right_Testis_mass) %>% pivot_longer(names_to = "Temp",  values_to = "Right_Testis", cols = c(Right_Testis1, Right_Testis_mass)) %>% 
  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Right_Testis) %>% distinct(Mouse_ID, .keep_all = T) 
SOTA <- full_join(SOTA, MT_Right_Testis) %>% select(-c(Right_Testis_mass, Right_Testis1))
rm(MT_Right_Testis)



## Tail_Length == "Tail_length", "LCd"
MT_Tail_Length <- SOTA %>% select(Mouse_ID, Tail_length, LCd) %>% pivot_longer(names_to = "Temp",  values_to = "Tail_Length", cols = c(Tail_length, LCd)) %>% 
  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Tail_Length) %>% distinct(Mouse_ID, .keep_all = T) 
SOTA <- full_join(SOTA, MT_Tail_Length) %>% select(-c(Tail_length, LCd))
rm(MT_Tail_Length)


## Trap_Date == "Capture"
MT_Trap_Date <- SOTA %>% select(Mouse_ID, Capture, Arrival) %>% pivot_longer(names_to = "Temp",  values_to = "Trap_Date", cols = c(Capture, Arrival)) %>% 
  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Trap_Date) %>% distinct(Mouse_ID, .keep_all = T) 
SOTA <- full_join(SOTA, MT_Trap_Date) %>% select(-c(Capture, Arrival))
rm(MT_Trap_Date)
## for 2015, there is only the Dissection date given -
## -> any idea how I can subtract 1 day from that date to make it conform with/and can include it in Trap_Date??


## Year == "Year", "year"
MT_Year <- SOTA %>% select(Mouse_ID, Year, year) %>% pivot_longer(names_to = "Temp",  values_to = "Year", cols = c(Year, year)) %>% 
  arrange(Mouse_ID) %>%  group_by(Mouse_ID) %>%  fill(c(everything()), .direction = "downup") %>%  ungroup() %>%  select(Mouse_ID, Year)  %>% distinct(Mouse_ID, .keep_all = T) 
SOTA <- full_join(SOTA, MT_Year) %>% select(-c(year)) %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 
rm(MT_Year)



#### MANUAL CORRECTIONS ########################################################
    ## correct Year for specific samples
SOTA$Year[ SOTA$Mouse_ID %in% c("SK_2903", "SK_2904")] <- 2014
SOTA$Year[ SOTA$Mouse_ID %in% c("AA_0330", "AA_0450", "AA_0451", "AA_0452")] <- 2017
    ## Status
SOTA <- SOTA %>% mutate(Status = replace(Status, Status == "(pregnant)", "pregnant"),
                                  Status = replace(Status, Status == "(young)", "young"))
    ## State
SOTA <- SOTA %>% mutate(State = replace(State, State == "Germany", "DE"),
                                  State = replace(State, State == "D", "DE"),
                                  State = replace(State, State == "Poland", "PL"))
    ## multiple Mice per Box (refers to the problem that faeces were not separated per individual)
SOTA <- SOTA %>% mutate(Multiple_Mice_per_Box = ifelse(Mouse_ID %in% c("AA_0514", "AA_0515", "AA_0513","AA_0349", "AA_0454", "ZZ_0037", "ZZ_0038"), TRUE, FALSE))

################################################################################
#### NARROW DOWN COLUMNS AND ROWS ##############################################

SOTA <- SOTA[colnames(SOTA) %in% c(basics, gen.loci, dissection.cols, oocyst.cols, initial.worms.cols)] %>%
  filter(!is.na(Longitude), 
         !is.na(Latitude))

################################################################################
#### 3. ADD NEW DATA ###########################################################

## new data and assays:
    ## new Dissection Data 
    ## Oocyst Counting (Flotation data)
    ## Infection Intensity (qPCR)
    ## Eimeria Genotyping
    ## Gene Expression (RT-qPCR)
    ## Immuno Data
        ## MES FACS
        ## CEWE Elisa (IFNy)
    ## Crypto (qPCR)
    ## data collected for non-mus rodents


#### 3.1. add new Dissection Data ####################################
Dis2018 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ18_Dissections.csv")
Dis2019 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ19_Dissections.csv")
Dis2021 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ21_Dissections.csv")

colnames(Dis2018)[colnames(Dis2018)%in%"Dissection_date"] <- "Dissection_Date"
colnames(Dis2018)[colnames(Dis2018)%in%"ASP"] <- "Aspiculuris"
colnames(Dis2018)[colnames(Dis2018)%in%"SYP"] <- "Syphacia"
colnames(Dis2018)[colnames(Dis2018)%in%"HET"] <- "Heterakis"
colnames(Dis2018)[colnames(Dis2018)%in%"MART"] <- "Taenia_martis"
colnames(Dis2018)[colnames(Dis2018)%in%"CP"] <- "Catenotaenia_pusilla"
colnames(Dis2018)[colnames(Dis2018)%in%"HM"] <- "Hymenolepis_microstoma"
colnames(Dis2018)[colnames(Dis2018)%in%"HD"] <- "Hymenolepis_diminuta"
colnames(Dis2018)[colnames(Dis2018)%in%"TM"] <- "Trichuris_muris"
colnames(Dis2018)[colnames(Dis2018)%in%"MM"] <- "Mastophorus_muris"
colnames(Dis2018)[colnames(Dis2018)%in%"Epididymis"] <- "Left_Epididymis"
colnames(Dis2018)[colnames(Dis2018)%in%"Body_weight"] <- "Body_Weight"
colnames(Dis2018)[colnames(Dis2018)%in%"Body_length"] <- "Body_Length"
colnames(Dis2018)[colnames(Dis2018)%in%"Tail_length"] <- "Tail_Length"
colnames(Dis2018)[colnames(Dis2018)%in%"Trap_date"] <- "Trap_Date"
colnames(Dis2018)[colnames(Dis2018)%in%"Feces"] <- "Feces_Weight"
colnames(Dis2018)[colnames(Dis2018)%in%"Embryo_left"] <- "Left_Embryo"
colnames(Dis2018)[colnames(Dis2018)%in%"Embryo_right"] <- "Right_Embryo"

colnames(Dis2019)[colnames(Dis2019)%in%"Dissection_date"] <- "Dissection_Date"
colnames(Dis2019)[colnames(Dis2019)%in%"ASP"] <- "Aspiculuris"
colnames(Dis2019)[colnames(Dis2019)%in%"SYP"] <- "Syphacia"
colnames(Dis2019)[colnames(Dis2019)%in%"HET"] <- "Heterakis"
colnames(Dis2019)[colnames(Dis2019)%in%"MART"] <- "Taenia_martis"
colnames(Dis2019)[colnames(Dis2019)%in%"CP"] <- "Catenotaenia_pusilla"
colnames(Dis2019)[colnames(Dis2019)%in%"HM"] <- "Hymenolepis_microstoma"
colnames(Dis2019)[colnames(Dis2019)%in%"HD"] <- "Hymenolepis_diminuta"
colnames(Dis2019)[colnames(Dis2019)%in%"TM"] <- "Trichuris_muris"
colnames(Dis2019)[colnames(Dis2019)%in%"MM"] <- "Mastophorus_muris"
colnames(Dis2019)[colnames(Dis2019)%in%"Epididymis"] <- "Left_Epididymis"
colnames(Dis2019)[colnames(Dis2019)%in%"Body_weight"] <- "Body_Weight"
colnames(Dis2019)[colnames(Dis2019)%in%"Body_length"] <- "Body_Length"
colnames(Dis2019)[colnames(Dis2019)%in%"Tail_length"] <- "Tail_Length"
colnames(Dis2019)[colnames(Dis2019)%in%"Trap_date"] <- "Trap_Date"
colnames(Dis2019)[colnames(Dis2019)%in%"Feces"] <- "Feces_Weight"
colnames(Dis2019)[colnames(Dis2019)%in%"Embryo_left"] <- "Left_Embryo"
colnames(Dis2019)[colnames(Dis2019)%in%"Embryo_right"] <- "Right_Embryo"

Dis2018 <- Dis2018[colnames(Dis2018) %in% c("Mouse_ID", basics, dissection.cols, initial.worms.cols, oocyst.cols)]
Dis2019 <- Dis2019[colnames(Dis2019) %in% c("Mouse_ID", basics, dissection.cols, initial.worms.cols, oocyst.cols)]
Dis2021 <- Dis2021[colnames(Dis2021) %in% c("Mouse_ID", basics, dissection.cols, initial.worms.cols, oocyst.cols)]

    ## merge
Dis1 <- full_join(Dis2018, Dis2019)
Dis2 <- full_join(Dis1, Dis2021)

SOTA$Left_Embryo <- as.integer(SOTA$Left_Embryo)
SOTA$Right_Embryo <- as.integer(SOTA$Right_Embryo)

SOTA <- full_join(SOTA, Dis2) %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 

    ## correct the Sex column
SOTA$Sex[grep("female*.", SOTA$Sex)] <- "F"
SOTA$Sex[grep("male*.", SOTA$Sex)] <- "M"


#### 3.2. add Oocyst Counting Data (Flotation Data) ############################
Eflot2018 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Eimeria_detection/HZ18_Eim_Flotation.csv")
Eflot2018$Ncells <- Eflot2018$Sume
Eflot2018$PBS_dil_in_mL <- Eflot2018$PBS_vol
Eflot2018$Feces_Weight <- Eflot2018$Feces
colnames(Eflot2018)[colnames(Eflot2018)%in%oocyst.cols]
Eflot2018 <- Eflot2018[colnames(Eflot2018)%in%c(basics,oocyst.cols)]

#### 3.3. add 2018 qPCR Data ###################################################
EqPCR2018 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Eimeria_detection/HZ18_qPCR.csv")
colnames(EqPCR2018)[colnames(EqPCR2018)%in%"delta"] <- "delta_ct_cewe_MminusE"
EqPCR2018 <- EqPCR2018[colnames(EqPCR2018)%in%c(basics, EqPCR.cols)]

#### 3.4. add 2018 Eimeria Genotyping Data #####################################
EimPCR <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Eimeria_detection/Svenja/table_ct_and_more.csv")
EimPCR$Mouse_ID <- gsub("CEWE_AA_", "AA_0", EimPCR$Name)
EimPCR$eimeriaSpecies <-  gsub("E\\. ", "E_", EimPCR$Eimeria.subspecies)
EimPCR$eimeriaSpecies[EimPCR$eimeriaSpecies%in%c("non infected", "Eimeria sp.")] <- "Negative"
EimPCR$Sex <- NULL
EimPCR <- EimPCR[colnames(EimPCR)%in%c(basics, EimGeno.cols)]

#### 3.5. add 2019 qPCR Data ###################################################
EqPCR2019 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Eimeria_detection/HZ19_CEWE_qPCR.csv")
colnames(EqPCR2019)[colnames(EqPCR2019)%in%"delta"] <- "delta_ct_cewe_MminusE"
colnames(EqPCR2019)[colnames(EqPCR2019)%in%"MC"] <- "MC.Eimeria"
EqPCR2019 <- EqPCR2019[colnames(EqPCR2019)%in%c(basics, EqPCR.cols)]


#### Merge
Detection18 <- merge(EimPCR, EqPCR2018)
Detection18 <- merge(Detection18, Eflot2018)
SOTA <- full_join(SOTA, Detection18) %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 
SOTA <- full_join(SOTA, EqPCR2019)   %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 


#### 3.6. add 2016-18 Gene Expression Data #####################################
Gene_Expression <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Gene_expression/HZ16-18_gene_expression.csv") %>% select(-c(X, HI)) 
Gene_Expression$Target <- gsub(pattern = "IL-6", replacement = "IL.6", x = Gene_Expression$Target)
colnames(Gene_Expression)[colnames(Gene_Expression)%in%"delta"] <- "delta_ct_cewe_MminusE"
colnames(Gene_Expression)[colnames(Gene_Expression)%in%"MC"] <- "MC.Eimeria"
Gene_Expression <- unique(Gene_Expression)
Gene_Expression <- Gene_Expression %>% pivot_wider(names_from = "Target", values_from = "NE")
SOTA <- full_join(SOTA, Gene_Expression)  %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 


#### 3.7 add 2019 CEWE Elisa ###################################################
CEWE_Elisa <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/HZ19_CEWE_ELISA.csv") %>% select(-X)
SOTA <- full_join(SOTA, CEWE_Elisa) %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 


#### 3.8 add 2019 MES FACS #####################################################
MES_FACS <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/HZ19_MES_FACS.csv") %>% select(-X)
SOTA <- full_join(SOTA, MES_FACS) %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 


#### 3.9 add 2019 Immuno #######################################################
Immuno19 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/HZ19_immuno.csv") %>% select(-X)
SOTA <- full_join(SOTA, Immuno19) %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T) 
colnames(Immuno19)[colnames(Immuno19)%in%"delta"] <- "delta_ct_cewe_MminusE"


#### 4. ADD CRYPTO DATA ########################################################
Crypto_qPCR <- read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Crypto_Detection.csv") %>% select(-X)
Crypto_qPCR <- Crypto_qPCR[colnames(Crypto_qPCR) %in% c(Crypto_qPCR.cols, "Mouse_ID")]
SOTA <- full_join(SOTA, Crypto_qPCR) %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T)


#### 5. ADD NON-MUS DATA #######################################################
Non_Mus <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Eimeria_detection/Other_rodents/rawdata_other_rodents.csv")
colnames(Non_Mus)[colnames(Non_Mus)%in%"Oocyst_g"] <- "OPG"
colnames(Non_Mus)[colnames(Non_Mus)%in%"Ocount_11"] <- "N_oocysts_sq1"
colnames(Non_Mus)[colnames(Non_Mus)%in%"Ocount_12"] <- "N_oocysts_sq2"
colnames(Non_Mus)[colnames(Non_Mus)%in%"Ocount_13"] <- "N_oocysts_sq3"
colnames(Non_Mus)[colnames(Non_Mus)%in%"Ocount_14"] <- "N_oocysts_sq4"
colnames(Non_Mus)[colnames(Non_Mus)%in%"Ocount_21"] <- "N_oocysts_sq5"
colnames(Non_Mus)[colnames(Non_Mus)%in%"Ocount_22"] <- "N_oocysts_sq6"
colnames(Non_Mus)[colnames(Non_Mus)%in%"Ocount_23"] <- "N_oocysts_sq7"
colnames(Non_Mus)[colnames(Non_Mus)%in%"Ocount_24"] <- "N_oocysts_sq8"
colnames(Non_Mus)[colnames(Non_Mus)%in%"Eimeria_ident"] <- "eimeriaSpecies"
colnames(Non_Mus)[colnames(Non_Mus)%in%"CEWE_Dct"] <- "delta_ct_cewe_MminusE"

    ## merge
SOTA <- full_join(SOTA, Non_Mus[colnames(Non_Mus) %in% c(basics, oocyst.cols, "delta_ct_cewe_MminusE", "eimeriaSpecies", EqPCR.cols, tissue.cols)])
SOTA <- SOTA %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T)


#### 6. add 2021 Dissection Data ###############################################
HZ21_Dis <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ21_Dissections.csv")
HZ21_Dis <- HZ21_Dis %>% mutate(Year = 2021)
Worms21 <- HZ21_Dis %>% select("Mouse_ID", 28:36)


  ## merge
SOTA <- full_join(SOTA, HZ21_Dis[colnames(HZ21_Dis) %in% c(basics, dissection.cols)])


  ## Non_Mus21 Data
Non_Mus21 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ21_Non_Mus.csv")
Non_Mus21 <- Non_Mus21 %>%
  mutate(Ticks = case_when(Ticks == T ~ T,
                           is.na(Ticks) ~ F),
         Fleas = case_when(Fleas == T ~ T,
                           is.na(Fleas) ~ F),
         Ectoparasites_Logical = case_when(Ectoparasites_Logical == T ~ T,
                                           is.na(Ectoparasites_Logical) ~ F),
         Year = 2021)

SOTA <- full_join(SOTA, Non_Mus21[colnames(Non_Mus21) %in% c(basics, dissection.cols, tissue.cols)])



#### MANUAL CORRECTION ######################################################
  # Worms
    ## Aspiculuris_sp
SOTA <- SOTA %>% mutate(Aspiculuris_sp =  case_when(!is.na(Aspiculuris_Syphacia) ~ Aspiculuris_Syphacia - Syphacia,
                                                          is.na(Aspiculuris_Syphacia) ~ Aspiculuris))

    ## Syphacia_sp
SOTA <- SOTA %>% mutate(Syphacia_sp = case_when(!is.na(Aspiculuris_Syphacia) ~ Aspiculuris_Syphacia - Aspiculuris_sp,
                                                is.na(Aspiculuris_Syphacia) ~ Syphacia))

    ## Trichuris == "Trichuris" "Trichuris_muris"
SOTA <- SOTA %>% select(-Trichuris)

    ## Taenia == Taenia taeniformis + Taenia martis
SOTA <- SOTA %>% 
  mutate(Taenia_sp = case_when( (Taenia == 0  & Taenia_martis == 0) | (Taenia == 0  & Taenia_taeniformis == 0) | Taenia == 0 | Taenia_martis == 0 ~ 0,
                                (Taenia == 1  & Taenia_martis == 1) | (Taenia == 1  & Taenia_taeniformis == 1) | Taenia == 1 ~ 1,
                                Taenia == 2 & Taenia_martis == 0 ~ 2,
                                Taenia == 3 & Taenia_martis == 3 ~ 3,
                                Taenia == 6 ~ 6,
                                Taenia_martis == 1 ~ 1,
                                Taenia_martis == 2 ~ 2,
                                Taenia_martis == 6 ~ 6))

    ## Heterakis
SOTA <- SOTA %>% select(-Heterakis_spumosa)
colnames(SOTA)[colnames(SOTA)%in%"Heterakis"] <- "Heterakis_sp"

    ## Hymenolepis
SOTA <- SOTA %>% mutate(Hymenolepis_sp = case_when(Hymenolepis_microstoma >= 0 ~ Hymenolepis_microstoma,
                                                   Hymenolepis_diminuta > 0 ~ Hymenolepis_diminuta))



SOTA <- full_join(SOTA, Worms21) %>% arrange(Mouse_ID) %>% group_by(Mouse_ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(Mouse_ID, .keep_all = T)


    ## Worms Presence
SOTA <- SOTA %>% mutate(Worms_presence = case_when(Aspiculuris_sp | Trichuris_muris | Taenia_sp | Heligmosomoides_polygurus | Heterakis_sp | Mastophorus_muris | Hymenolepis_sp | Catenotaenia_pusilla > 0 ~ T))

 
#### 7. SELECT NEEDED COLUMNS ##################################################
SOTA <- SOTA[colnames(SOTA) %in% c(basics,
                                   tissue.cols,
                                   Crypto_qPCR.cols,
                                   dissection.cols,
                                   EimGeno.cols,
                                   EqPCR.cols,
                                   gen.loci,
                                   Gene.Exp.cols,
                                   oocyst.cols,
                                   #initial.worms.cols,
                                   final.worms.cols)]
                                            
SOTA[colnames(SOTA) %in% c(basics, final.worms.cols)] %>% group_by(Year) %>% vis_miss()


write.csv(SOTA, "SOTA_Data_Product.csv")



