# Eimeria_Wild, Data analysis for parasites in the house mouse hybrid zone  

###  Current Maintainer: Lubomir Bednar

# <span style="color:darkred">1. Organisation<span>

Mouse_Eimeria_Datababasing


This repository is for storage of all !relevant material! on wild mice capture and Eimeria detection. Each file should be named according to the tamplate of:
HZyear_Parasite_DetectionMethod_OtherInfo.format (Detection assays)
HZYear_Trap (Trapping information)
HZYear_Dissection (Dissection information)
HZYear_Genotypes (Genotype information on captured mice)
The templates for filling each mandatory column is indicated next to the column names.

E.g.: HZ18_Eimeria_qPCR_raw.csv
E.g.: HZ18_Trap
E.g.: HZ18_Dissection
E.g.: HZ18_Genotypes

Folders:
data = Main folder containing Cryptosporidium, Eimeria_detection, Field_data
	
	Cryptosporidium =  A folder dedicated to detection of Cryptosporidium parasites, intended for storage of !raw! data. Must contain in the very least a "Mouse_ID"(AA_xxxx) and "Year"(xxxx) columns.

	Eimeria_detection = Like Cryptosporidium but for Eimeria. All detection methods and flotation raw data go in here to be used by the group. Must contain in the very least a "Mouse_ID" and "Year" columns.

	Field_data = This is the field trip data repository for Trap tables, Dissection tables and Genotype tables. In addition a Merge table can be made combining all 3 to link our data from one year together 		in a complete table. Each type of table has it's own prerequisites (Merge inherits prerequisites from it's components). 
	
	Trap tables must include: "Latitude"(xx.xxxx), "Longitude"(xx.xxxx), "Trap_set_date" (dd/mm/yy), "Trap_collect_date(dd/mm/yy), "Trap_set_time"(hh:mm), "Trap_collect_time"(hh:mm) "Address", , 		"Trap_number", "Team", "Mus_caught", "Rodents_caught" and "Comments".
	
	Dissection tables must include: "Mouse_ID" (AA_xxxx), "Address", "Latitude"(xx.xxxx), "Longitude"(xx.xxxx), "Trap_date"(dd/mm/yy), "Dissection_date"(dd/mm/yy), "Sex", "Status", "Species"(Genus_species), 		"Ectoparasites" (Logical), "Body_weight"(xx.xx), "Body_length"(xx.x), "Tail_length"(xx.x), "Spleen"(x.xxxx), "Left_testis"(x.xxxx), "Right_testis"(x.xxxx), "Epididymis"(x.xxxx), 	  	 		"Seminal_vesicle"(x.xxxx), "Embryo_left"(quantity), "Embryo_right"(quantity), "Feces"(x.xx), "ASP", "SYP", "HET", "MART", "CP", "HD", "HM", "MM" and "TM" (Worm acronyms, numeric).
	
	Genotype tables mus include: "Mouse_ID"(AA_xxxx, "Sex", "Latitude"(xx.xxxx), "Longitude"(xx.xxxx) and "HI"(x.xx+).


R = This folder contains data processing scripts for cleaning the raw data tables and analyzing them once cleaned. Here you can find all the necessary code for how previous work was processed and help yourself to useful functions to make your analysis easier and !replicable!. The case is same here as in Eimeria_Lab = Please write your script in a way to load the data from raw. GitHub files. Name your scripts by the template + whatever they are made for.

	e.g.: 
	#load libraries for loading raw github files
	library(httr)
	library(RCurl)
	#read in as URL from GitHub
	cell.countsURL <- "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/3_recordingTables/E7_112018_Eim_FACS_cell_counts_processed.csv"
	cell.counts <- read.csv(text = getURL(cell.countsURL)) 
	................................................................
	#after processing and cleaning, write to the appropriate repo and Git it
	write.csv(E7, "./Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_FACS_clean.csv", quote = FALSE)

.git = A folder for GitHub use when initializing new repo or cloning an existing one. I wouldn't touch it unless you know what you're doing.

DO NOT edit the raw data in these folders. You can clone the repos, work with the data and generate new tables as you please. The envisioned structure is that each manuscript has it's own folder/repo, the raw data is loaded from the raw. GitHub files and saved in the manuscript folder. Then all can be edited and analysed there.
 "makeMiceTable.R" should be updated. Never touch the raw data, they are our *precioussss*
