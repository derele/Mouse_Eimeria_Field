# Storage and curation of data generated from the European house mouse hybrid zone

This repository is for storage of all relevant data on wild mice
capture and 

including detection of Eimeria, Cryptosporidium, worms,
detection methods and gene expression.

# Structure:

## [data_access_code](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data/data_access_code) = R code for growing the combined final data products

### [Field_data](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data/Field_data) = Trapping, Dissections and Genotype data. The first two taken during the field trip the latter (genotyping) from Jaroslav Pialek.

### [Eimeria_detection](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data/Eimeria_detection) = Eimeria qPCRs and flotation result

For documentation of the assays behind datasets and the standardised column names see !ADD LINK!

### [Gene_exprssion](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data/Gene_expression) = Gene expression data

For documentation of the assays behind datasets and the standardised column names see !ADD LINK!

### [Cryptosporidium](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data/Cryptosporidium) = Cyptosporidium screening

For documentation of the assays behind datasets and the standardised column names see !ADD LINK!

## [data_access_code](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data/data_access_code) = R code for growing the combined final data products

For the current of the art see !ADD LINK!

## [data_creation_code](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data/data_creation_code) = R for the processing of raw data before storage in [data](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data).

## [data](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data) = contains all data from field excursions and follow up lab work it also has the templates for the table that should be use for recording 


# 1. Accessing data:
## 1.1. General description:

All data in this repository has been processed and saved as a clean
table according to the corresponding
[template](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data/Templates)

All file names contain information to distinguish the year of
collection, assay (where applicable), type of fieldwork (where
applicable) and a format.  In addition, elements of the name may
designate: [tissue
type](https://github.com/derele/Mouse_Eimeria_Field/tree/master/Tissue_labels.csv)
or assay type (e.g., RT-qPCR)

Exmaple: HZ19_CEWE_Eim_qPCR.csv This means the table contains
information generated from 2019, the tissue used in the assay was
Caecum, assay was targeted to screen for Eimeria, the assay type was a
qPCR and the table is in a .csv format.

## 1.2. Examples:
### 1.2.1. Example 1 

### 1.2.1. Example 2 

# 2. Adding data:
## 2.1. General description:
Each file should be named according to the tamplate of:
HZYear_TypeOfFieldwork.format
E.g.: HZ19_Dissections.csv
or
HZYear_TissueType_Parasite_AssayType.format
E.g.: HZ19_CEWE_Eim_RT-qPCR.csv

Raw data should be stored here and processed using code saved here as well. Both should be subsequently deleted once a clean table exists. The raw data and code should be both commited and pushed to git to keep track of events. Commit messages should contain information on what files are being handled.

General rule is:
1. upload raw data table
2. upload code to process raw data table
3. upload clean data table
4. delete raw data table and code

## 2.2. Examples:
### 2.2.1. Adding genotype data

### 2.2.2. Adding qPCR data
