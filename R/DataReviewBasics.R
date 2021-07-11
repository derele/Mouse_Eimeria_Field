### We start with the structure of the file
### "data/MiceTable_fullEimeriaInfos_2014to2017.csv"

ALL <- read.csv("data/MiceTable_fullEimeriaInfos_2014to2017.csv")

colnames(ALL)
### this has 141 columns! Crazy! The way Alice had organized this is
### to leave every column in the data, even the ones with different
### spelling (produced by students or collaboration partners) are
### still in the data! So all data in this table is for sure  still there!

###  We now need work with this and need to be selective. Here are
###  some selections of the main columns that could be helpful. For
###  other tasks you might be other columns helpful.  (for the ones
###  with abberant spelling ("Year/year", "Latitude/latitude", etc..)
###  I checked that we have the "main column" with all/most of the
###  information

basics <- c("Mouse_ID", "Sex", "Longitude", "Latitude", "Year")

gen.loci <- c("mtBamH", "YNPAR", "X332", "X347", "X65", "Tsx", "Btk", "Syap1",
              "Es1", "Gpd1", "Idh1", "Mpi", "Np", "Sod1", "Es1C", "Gpd1C",
              "Idh1C", "MpiC", "NpC", "Sod1C", "Zfy2", "SRY1", "Y", "HI_NLoci",
              "HI")

parasite.cols <- c("Aspiculuris_tetraptera", "Syphacia_obvelata", "Trichuris_muris",
                   "Taenia_taeniformis", "Flea", "Mix_Syphacia_Aspiculuris",
                   "Heterakis_spumosa", "Mastophorus_muris", "Hymenolepis_microstoma",
                   "Catenotaenia_pusilla", "Cysticercus", "Ectoparasites",
                   "Worms_presence", "Hymenolepis_diminiuta", "Taenia_martis",
                   "Heligmosomoides_polygurus", "Taenia", "Aspiculuris_Syphacia",
                   "Trichuris", "Heterakis", "Mastophorus")

oocyst.cols <- c("counter", "Feces_g", "Date_count", "N_oocysts_sq1",
                 "N_oocysts_sq2", "N_oocysts_sq3",  "N_oocysts_sq4",
                 "N_oocysts_sq5", "N_oocysts_sq6", "N_oocysts_sq7",
                 "N_oocysts_sq8", "mean_neubauer", "PBS_dil_in_mL", 
                 "OPG", "Ncells")

EqPCR.cols <- c("delta_ct_ilwe_MminusE", "delta_ct_cewe_MminusE",
                ## from 2018 on AN IMPORTANT IMPROVEMENT!!!
                "MC.Eimeria")

EimGeno.cols <- c("n18S_Seq", "COI_Seq", "ORF470_Seq", "eimeriaSpecies")


### Now we add the main 2018 data! #######################################

################### 2018 #################################################

### The basics from the dissection
Dis2018 <- read.csv("data/Field_data/HZ18_Dissections.csv")

### The mouse genotyping from Jarda's table
Gen2018 <- read.csv("data/Field_data/HZ18_Genotypes.csv")

## So Jarda's genotyping table acutally has all the main data!  As
## long as we don't want to look at the other parasites we have to do
## nothing to this befor the merge
DisGen2018 <- Gen2018[, colnames(Gen2018)%in%c(basics, gen.loci), ]

### So adding the oocysts, the qPCR and the Eimeria Species genotyping...
Eflot2018 <- read.csv("data/Eimeria_detection/HZ18_Eim_Flotation.csv")

colnames(Eflot2018)
## Flotation columns have different names, check them if you need and
## adjust (translate them into proper __ oocyst.cols __ ... for now
## we're believing just in the correctness of the OPG
## calculation... and add only OPG (which is the only proper __
## oocyst.cols __

### for now simply Eflot2018[, colnames(Eflot2018)%in%oocyst.cols]
### will add only the OPG: see
colnames(Eflot2018)[colnames(Eflot2018)%in%oocyst.cols]

EqPCR2018 <- read.csv("data/Eimeria_detection/HZ18_qPCR.csv")
## the "delta" is __ delta_ct_cewe_MminusE __ as we only do cecum wall
## (we don't do illeum as vermiformis is so rare).
colnames(EqPCR2018)[colnames(EqPCR2018)%in%"delta"] <- "delta_ct_cewe_MminusE"

## now we will add
colnames(EqPCR2018)[colnames(EqPCR2018)%in%EqPCR.cols]

