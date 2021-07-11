# Storage and curation of data generated from the European house mouse hybrid zone

This repository is for storage, curation and documentation of all
relevant data on wild mice capture follwing laboratory assays
including detection of Eimeria, Cryptosporidium, worms, gene
expression, etc...

# Structure:

## [data_products](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data_products) 

Analysis ready datasets to be used (or which have bin used) in
research projects:

- [MiceTableMusAliceArticle.csv](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data_products/MiceTableMusAliceArticle.csv): Has been compiled using the script [MiceTableMusAliceArticle.R](https://github.com/derele/Mouse_Eimeria_Field/blob/master/R/input2product/MiceTableMusAliceArticle.R)


## [data_input](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data_input) 

###  [Mouse_data](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data_input/Mouse_data) 

This contains two important datasets for each year: 

- Dissections 

Trapping data is originally compiled in as "trapping data"
(HZ\\d\\d_Trap.csv) by the catching teams. We also store this
farm-level data here to be able to resolve potential potential
problems in the data later. The trapping data, however, is already
merged with dissection data during the field trip to control the
accuracy of locations and assignment of mice to those on the
spot. This results in combined mouse dissection datasets, which follow
the the nameing scheme : HZ(\\d\\d for year)_Dissections.csv,
e.g. HZ18_Dissections.csv

For documentation on the trapping procedure and dissection (including
standardised column name) see: !ADD LINK!


- Genotypes 
Genotype data is compiled for us by Jaroslav Pialek (for
2014-2019). It follows the same naming scheme: e.g. HZ18_Genotypes.csv

For documentation on the assays behind the genotyping datases
(including standardised column name) see:: !ADD LINK!


### [Eimeria_detection](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data_input/Eimeria_detection) 

We screen our samples for parasite infections during the following
labowork, this results in two datasets Eimeria qPCRs and flotation result

For documentation of the assays behind datasets and the standardised
column names see !ADD LINK!

WHAT about EIMERIA species TYPING!? This data seems missing or not
added in a standardized way!


### [Gene_expression](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data_input/Gene_expression) 

Datasets are: 
- x
- Y
- Z

For documentation of the assays behind datasets and the standardised
column names see !ADD LINK!

### [Cryptosporidium](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data_input/Cryptosporidium)

For documentation of the assays behind datasets and the standardised column names see !ADD LINK!


### Templates:

[template](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data_input/Templates)

All file names contain information to distinguish the year of
collection, assay (where applicable), type of fieldwork (where
applicable) and a format.  In addition, elements of the name may
designate: [tissue
type](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data_input/Templates/Tissue_labels.csv)
or assay type (e.g., RT-qPCR)

Exmaple: HZ19_CEWE_Eim_qPCR.csv This means the table contains
information generated from 2019, the tissue used in the assay was
Caecum, assay was targeted to screen for Eimeria, the assay type was a
qPCR and the table is in a .csv format.

## [R](https://github.com/derele/Mouse_Eimeria_Field/tree/master/R) 

Code required for data processing and curation.

### [input2product](https://github.com/derele/Mouse_Eimeria_Field/tree/master/R/input2product)

Code producing final data producsts in
[data_product](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data_product),
currently:

- [DataReviewBasics.R](https://github.com/derele/Mouse_Eimeria_Field/tree/master/R/input2product/DataReviewBasics.R):

 A comprehensive review of "all*" the data throughout the years. From
 a huge candidate list of potential columns, this estabishes some
 groups of standardized meaningful columns.  Pre-2017 data is based on
 [MiceTableMusAliceArticle.csv](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data_products/MiceTableMusAliceArticle.csv)
 in which Alice had compiled the most coprehensive list of possible
 columns. *We select a workable subset of columns (as documented for
 individual datasets, see abofe) and add 2018 and 2019 data.
 
 We will expand the workable subset of columns further...  
 

- [MiceTableMusAliceArticle.R](https://github.com/derele/Mouse_Eimeria_Field/tree/master/R/input2product/MiceTableMusAliceArticle.R):
  this has been used to compile
  [MiceTableMusAliceArticle.csv](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data_products/MiceTableMusAliceArticle.csv).It's
  currently not executable in the present re-structed
  repository. (Also, as a coding advice: please don't write the whole
  content of a scritp into one single function)


### [raw2input](https://github.com/derele/Mouse_Eimeria_Field/tree/master/R/raw2input) 

Contains code that might be helpful to curate the data _before_ it
goes into
[data_input](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data_input),
currently: 

- [HMHZ_Functions.R](https://github.com/derele/Mouse_Eimeria_Field/tree/master/R/raw2input/HMHZ_Functions.R):
  a few functions to e.g. translate degree into decimal lat/long and
  some ideas on clustering corrdinates for farm-level analysis (might
  be rather analysis code).


### [analysis](https://github.com/derele/Mouse_Eimeria_Field/tree/master/R/analysis) 

Contains some old data analysis scipts. PLEASE DON'T use this
repository TO TRACK DATA ANALYSIS projects. The data analysis scripts
tracke here should be purely for reviewing the correctnes of data. The
older scirpts are left here for inspiration for now. We might remove
them. 





### Luke's documentation of his procedure read

Raw data should be stored here and processed using code saved here as
well. Both should be subsequently deleted once a clean table
exists. The raw data and code should be both commited and pushed to
git to keep track of events. Commit messages should contain
information on what files are being handled.

General rule is:
1. upload raw data table
2. upload code to process raw data table
3. upload clean data table
4. delete raw data table and code

We need to figure out what this means for each individual (still
undocumented dataset)
