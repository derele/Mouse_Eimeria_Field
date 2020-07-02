# Mouse_Eimeria_Field, storage of data generated from the European house mouse hybrid zone (EHMHZ)  
This repository is for storage of all relevant data on wild mice capture, including detection of Eimeria, Cryptosporidium, worms, detection methods and gene expression. 

# Structure:
## [data](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data) = contains all cleaned up data generated during our field excursions and 
templates for tables

### [Field_data](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data/Field_data) = records from the field excursions, contains Trapping,
Dissections and Genotype data.

### [Eimeria_detection](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data/Eimeria_detection) = clean tables of Mouse-Eimeria qPCRs and
flotation result tables.

### [Templates](https://github.com/derele/Eimeria_Lab/tree/master/data/Templates) = examples of what correspondingly named tables should look like

#### [Gene_exprssion](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data/Gene_expression) = Folder for storing tidy gene expression 
tables from the wild

### [Cryptosporidium](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data/Cryptosporidium) = Similar to Eimeria_detection. For storage 
of tidy tables related to Cyptosporidium screening in the wild.

## [data_access_code](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data/data_access_code) = contains examples of R code related to accessing and reading information from the [data](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data) folder.

## [data_creation_code](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data/data_creation_code) = contains examples of R code related to processing raw data and making it suitable to be in the [data](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data) folder.


# 1. Accessing data:
## 1.1. General description:

All data in this repository has been processed and saved as a clean table according to the corresponding 
[template](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data/Templates)

All file names contain information to distinguish the year of collection, assay (where applicable), type of fieldwork (where applicable) and a format.
In addition, elements of the name may designate: [tissue type](https://github.com/derele/Mouse_Eimeria_Field/tree/master/Tissue_labels.csv) or assay type (e.g., RT-qPCR)

Exmaple: HZ19_CEWE_Eim_qPCR.csv
This means the table contains information generated from 2019, the tissue used in the assay was Caecum, assay was targeted to screen for Eimeria, 
the assay type was a qPCR and the table is in a .csv format.

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

Raw data should be stored here and processed using code saved here as well. Both should be subsequently deleted once a clean table exists.
The raw data and code should be both commited and pushed to git to keep track of events. Commit messages should contain information on what
files are being handled.

General rule is:
1. upload raw data table
2. upload code to process raw data table
3. upload clean data table
4. delete raw data table and code

## 2.2. Examples:
### 2.2.1. Adding genotype data

### 2.2.2. Adding qPCR data
