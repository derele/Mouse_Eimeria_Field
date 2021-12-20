---
title: "README-2.0"
output:
  html_document:
    toc: yes
    toc_float: yes
    number_section: yes
    fig_caption: true
  pdf_document:
    toc: yes
  always_allow_html: true
date: "`r format(Sys.time(), '%d %B %Y')`"
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(fig.align = 'center', echo = TRUE)
```

``` {r libraries-and-data, echo = F, message = F, warning = F}
library(tidyverse)
library(knitr)
```

# Storage and Curation of data generated from the HMHZ

This repository is for storage, curation and documentation of all
relevant data on wild mice capture following laboratory assays
including detection of _Eimeria_, _Cryptosporidium_, worms, gene
expression, etc...



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
Main Documentation file in csv format. New variables should be added continuously.

``` {r Documentation-MEF, echo = F, message = F, warning = F, include = TRUE} 

Docu_MEF <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/documentation/Documentation_MEF1.csv")
kable(Docu_MEF[,1:7])

```