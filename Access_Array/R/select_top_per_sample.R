## get the proper location for the functions (mainly readEseq)
source("/data/Parascaris/scripts/functions.R")

files <- list.files("/data/RAW_SEQ/2015_16_2015_18/qual_filter/sorted_amps_trimmed/flashed/" , "*.fasta")
setwd("/data/RAW_SEQ/2015_16_2015_18/qual_filter/sorted_amps_trimmed/flashed/")

top.samples <- list()

for (f in files) { 
    seq <- read.ESeq(f)
    seq <- seq[!grepl("foo", names(seq))]
    samples <- gsub("RECORD:(.*?_chip\\d)_.*", "\\1", names(seq))
    top <- by(cbind(samples, seq), samples, function (x) {
                  seq.table <- table(x$seq)
                  hack <- names(seq.table[seq.table == max(seq.table)])[[1]]
                  hack
              })
    tip.top <- rbind(unlist(top))
    top.samples[[f]] <- tip.top
}

ALL <- melt(top.samples)

##  in one file for  blast screeening

all.seq.one <- ALL$value
loc.names <- strsplit(ALL$L1, "_")
loc.names <- unlist(lapply(loc.names, function(x){
                               paste(x[c(1, 2)], collapse="_")
                           }))
names(all.seq.one) <- paste(loc.names, ALL$X2, sep=";")

## change only when improving the method

## write.ESeq(all.seq.one, "/data/RAW_SEQ/2015_16_2015_18/qual_filter/sorted_amps_trimmed/flashed/top_per_sample/blastscreen/all_named.fasta")

## do the blasts and balst2alltax
vs_ef <- read.csv("/data/RAW_SEQ/2015_16_2015_18/qual_filter/sorted_amps_trimmed/flashed/top_per_sample/blastscreen/all_named_vs_Efal.taxtab", as.is=TRUE)
vs_ef$species <- "Eimeria falciformis"
vs_ef$genus <- "Eimeria"
vs_ef$family <- "Eimeriidae"
vs_ef$phylum <- "Apicomplexa"
vs_ef$kingdom <- "undef"
vs_ef$superkingdom <- "Eukaryota"

vs_nt <- read.csv("/data/RAW_SEQ/2015_16_2015_18/qual_filter/sorted_amps_trimmed/flashed/top_per_sample/blastscreen/all_named_vs_nt.taxtab",
                  as.is=TRUE)

## nr seems not necessary
## vs_nr <- read.csv("/data/RAW_SEQ/2015_16_2015_18/qual_filter/sorted_amps_trimmed/flashed/top_per_sample/blastscreen/all_named_vs_nr.taxtab")

vs_tax<- rbind(vs_ef, vs_nt)
vs_tax <- vs_tax[order(vs_tax$query, vs_tax$bitscore),]

vs_tax <- vs_tax[!duplicated(vs_tax$query),]
exclude <- vs_tax$query[vs_tax$phylum!="Apicomplexa"]

Eimeria.seq.one <- all.seq.one[!names(all.seq.one)%in%exclude]
Eimeria.seq.one <- Eimeria.seq.one[!grepl("H20", names(Eimeria.seq.one))]
Eimeria.seq.one <- Eimeria.seq.one[!grepl("H2O", names(Eimeria.seq.one))]
Eimeria.seq.one <- Eimeria.seq.one[!grepl("Agave", names(Eimeria.seq.one))]

split.names <- strsplit(names(Eimeria.seq.one), ";")

loci <- unlist(lapply(split.names, "[", 1))
Eimeria.seq.one <- as.character(Eimeria.seq.one)

names(Eimeria.seq.one) <- unlist(lapply(split.names, "[", 2))

## by(cbind(loci, seq=Eimeria.seq.one), loci,  function (x){
##     filename <- paste("/data/RAW_SEQ/2015_16_2015_18/qual_filter/sorted_amps_trimmed/flashed/top_per_sample/blastscreen/",
##                       x$loci[[1]], "_screened.fasta", sep="")
##     write.ESeq(x$seq, filename)}
##    )
    
## run muscle outside R
## parallel muscle -clwstrict -in {} -out aln/{.}.aln ::: *.






