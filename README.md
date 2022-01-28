# Storage and Curation of data generated from the HMHZ

This repository is for storage, curation and documentation of all
relevant data on wild mice capture following laboratory assays
including detection of _Eimeria_, _Cryptosporidium_, worms, gene
expression, etc...

# Structure:

## [data_input](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data_input)
All input data is kept together (can be merged by the identifier)
_"Mouse_ID"_, other columns should not be repeated in input data
tables. The columns are documented for the different categories of input 
data as indicated in the documentation folder.


## [data_products](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data_products)
Analysis ready datasets to be used (or which have been used) in
research projects


## [R](https://github.com/derele/Mouse_Eimeria_Field/tree/master/R)
Code required for data processing and curation

### [raw2input](https://github.com/derele/Mouse_Eimeria_Field/tree/master/R/raw2input)
  Contains code that might be helpful to curate the data _before_ it
  goes into
  [data_input](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data_input)

### [input2product](https://github.com/derele/Mouse_Eimeria_Field/tree/master/R/input2product)
  Code producing final data products in [data_products](https://github.com/derele/Mouse_Eimeria_Field/tree/master/data_products)

### [analysis](https://github.com/derele/Mouse_Eimeria_Field/tree/master/R/analysis)
Contains some old data analysis scripts. PLEASE DON'T use this
repository TO TRACK DATA ANALYSIS projects. The data analysis scripts
tracks here should be purely for reviewing the correctness of data. The
older scripts are left here for inspiration for now. We might remove
them.

## [documentation](https://github.com/derele/Mouse_Eimeria_Field/tree/master/documentation)
[Main Documentation file](https://github.com/derele/Mouse_Eimeria_Field/blob/master/documentation/Documentation_MEF.csv) in csv format. New variables should be added continuously.

- Assay: What kind of assay the variable is derived from
- Variable: Name in data tables
- Description: precise description of variable meaning
- Type: Data type (numeric, logical, character)
- Unit: if applicable
- Values: for orientation
- Correction needed: specifies that data in [SOTA](https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/SOTA_Data_Product.csv) may need some modification because of comma-errors or specific issues (date unification, different input across the years)
- aka: outdated names

## [Freezer Samples (-80)]
[Here](https://docs.google.com/spreadsheets/d/1AqSXkeK1bVVAKrHtJ3BtzcNf---XHzShCdP6IQAZI_Y/edit?usp=sharing) you can find the link to all freezer samples. If you take boxes upstairs to the -20 freezer, please adjust the sheet accordingly, so we are always up to date!
