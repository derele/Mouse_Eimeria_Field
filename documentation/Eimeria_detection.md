## Oocyst flotation

Add here... 

We report these in HZ\\d\\d_qPCR.csv, standard columns are: 

- counter: character, Person who counted
- Feces_g: numeric, amount of feces used in flotation
- Date_count: date counted
- N_oocysts_sq1 ... sq8: numeric, individual count for each single
   square on the neubaure chamber (up to 8)
- mean_neubauer: numeric, mean of the 8 squares 
- PBS_dil_in_mL: numeric, in which volume of PBS whas
- Ncells: number of neubauer "cells" (squares) counted 
- OPG: oocysts per gram feces (calculated from the above)

In some cases only OPG data might be available or rather data would
need to be re-formatted to access all the raw values.

## Eimeria detection qPCR

We perform relative qPCR for detection and quantification of Eimeria
DNA in intestinal tissue. We amplify a locus in the nuclear genome of
the house mouse and a locus in the mitochondrial genome (COI) of
Eimeria. We then calculate a "delta" between the two ct values.


We report these in HZ\\d\\d_qPCR.csv, standard columns are: 

- delta_ct_ilwe_MminusE: threshold cycle for mouse minus Eimeria in
  Ileum tissue. Only E. vermiformis is (at low pervalence in Ileum
  tissue) and we therfore don't obtain this data-type for all years.  

- delta_ct_cewe_MminusE: threshold cycle for mouse minus Eimeria in
  Caecume tissue. E. ferrisi and E. falcifromis are detected here. We
  should have this (as coprehensively as possible) for every year!

- MC.Eimeria: TRUE/FALSE. This was established in 2018 as an
  improvement over the '> -5 delta ct rule' for identification of
  Eimeria -positive qPCRs. Melting curves have to show a drastic drop
  at XX°C to indicate melting of a proper Eimeria COI amplification
  product. It might be added where possible for per 2018 data post-hoc
  (if melting curves exist for a review of raw data). 


Files with non standard formatting are currently:

- HZ18_qPCR.csv (CEWE data)
- HZ19_CEWE_qPCR.csv (CEWE data with additional mouse and Eimeria CTs)
