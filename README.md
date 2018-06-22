# Eimeria_Wild, Data analysis for parasites in the house mouse hybrid zone  

###  Maintainers: Alice Balard, Emanuel Heitlinger
### last update : 22 june 2018

# <span style="color:darkred">1. Organisation<span>

* figures
* R : scripts in the main folder, *functions* in a sub folder. Put there only functions that should be called by scripts.
* raw_data : data on Eimeria, regularly updated (flotation, PCR, qPCR)

# <span style="color:darkred">2. Important things to know<span>

The mice data from the HMHZ should not be touched. They are in our team Seafile folde "Data_important". You can put it in your local computer, and to generate the mice table that contains all the data, you should do the following:

> miceTable <- makeMiceTable("path/to/Data_important/")

Each year, when new data are generated (HI, etc.), the R script "makeMiceTable.R" should be updated. Never touch the raw data, they are our *precioussss*
