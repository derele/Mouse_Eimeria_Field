## We only use BioConductor library, which we load beforehand
library(rtracklayer)
## and R's most fancy plotting library
library(ggplot2)

## all other functions in this code are hand crafted and you can check
## them in the file functions.R We execute this file to make the
## functions available to us
source("functions.R")

falGenome <- read.ESeq("Efal_final.fa")
verGenome <- read.ESeq("Ever_final.fa")
nieGenome <- read.ESeq("Enie_final.fa")

## lets see how big our genome assemblies are:
sum(nchar(falGenome))

## Now we read the annotation data 
gff <- import.gff3("Efal_genes.gff3")
## From this we get only coding sequence exons
gff.cds <- subset(gff, gff$type%in%"CDS")
## Can you describe what you see in that data?

## We export it to a regular data frame to work with basic R
gff.cds.df <- as.data.frame(gff.cds)

## We extract the corresponding sequence from the Eimeria genome
falExome <- unlist(lapply(1:nrow(gff.cds.df), function(i){
    substr(falGenome[as.character(gff.cds.df[i, "seqnames"])],
           gff.cds.df[i, "start"],
           gff.cds.df[i, "end"])
}))

## and give it ids (names) derived from its positions on genomics
## contigs
names(falExome) <- paste(gff.cds.df$seqnames, ":" ,
                         gff.cds.df$start, "-",
                         gff.cds.df$end, sep = "")

## How bit is that exome? (Exome = all Exons)
sum(nchar(falExome))

## But we want to amplify also from "closely" related species. We will
## use BLAST to compare. Let's use only those exons which
## are larger than 120bases.
## write.ESeq(falExome[nchar(falExome)>120], "Efal_exons120.fasta")

## BLAST searches have been executed for you in the shell (comand
## line). The commands looked like this:

## blastn -query Efal_exons120.fasta -db
## ../nieschulzi/nieschulzi_final.fa -outfmt 6 -evalue 1e-5
## -max_target_seqs 1 -num_threads 24 > Fal_vs_nie.blt

##  blastn -query Efal_exons120.fasta -db
##  ../vermiformis/vermiformis_final.fa -outfmt 6 -evalue 1e-5
##  -max_target_seqs 1 -num_threads 24 > Fal_vs_ver.blt

## We now read the blast reports (check the function to do so in the
## functions.R file).
falVSver <- read.blast6("Fal_vs_ver.blt")

## rows of the table
nrow(falVSver)
## the predicted exons fall on that many E. vermiformis contigs:
length(unique(falVSver$subject_id))

### look how many matches are longer than 120bp and over 80% idendity
summary.factor(falVSver$alignment_length>120)
summary.factor(falVSver$identity>80)
## by far most of them!

### select only those matches longer than 120 and over 80% idendity
falVSver <- falVSver[falVSver$alignment_length>120 &
                     falVSver$identity>80 , ]

falVSnie <- read.blast6("Fal_vs_nie.blt")
### select matches longer than 120bp and over 80% idendity,
### alternative way to write it
falVSnie <- subset(falVSnie,
                   falVSnie$alignment_length>120 &
                   falVSnie$identity>80 )

## lets see how the identity values are distributed
id4plot <- rbind(cbind(falVSnie, species="E.nieschulzi"),
                 cbind(falVSver, species="E.vermiformis"))
  
png("identities.png", width=4, height=4, units='in', res=300)
ggplot(id4plot, aes(identity, color=species)) + geom_density() +
  scale_x_continuous("identity in %") +
  scale_color_discrete(guide = guide_legend()) +
  theme(legend.position="bottom")
dev.off()


## select the sequence parts of the vermiformis and nieschulzi genome,
## which are similar to the Efalciformis exome
verExome <- get.Start.End.blast(falVSver, verGenome)

## same "Exome selection" for the other genome
nieExome <- get.Start.End.blast(falVSnie, nieGenome)

## Write the "Exomes" to fasta files:
## write.ESeq(verExome, "verExome.fasta")
## write.ESeq(nieExome, "nieExome.fasta")

## I have performed a all-vs-all BLAST of the three exomes for you and
## clustered the output with MCL http://micans.org/mcl/.  We read the
## output into a list like this
mcl <- readLines("out.Eeexon.mci.I40")
mcl <- lapply(mcl, strsplit, "\t")
mcl <- lapply(mcl, function(x) x[[1]])

## We already exclude list elements with less than 3 entries at this
## point. They can't have entries from all three species:
mcl <- mcl[unlist(lapply(mcl, length))>2]

## We tabulate the number of sequences from each species in each
## cluster
Tmcl <- tabulate.mcl(mcl)

## Clusters should contain exons from all three species for but not
## more than 1 from each
select.cluster <- rowSums(Tmcl)==3

## If a species had no member in a cluster we got NA in the table and
## get will use FALSE now in the selection to not select the cluster
select.cluster[is.na(select.cluster)] <- FALSE

## This tells you how many clusters we will discarded:
summary.factor(select.cluster)

## Tis is the top and bottom of the selected table listing the numbers
## for species per cluster
head(Tmcl[select.cluster,])
tail(Tmcl[select.cluster,])

Stmcl <- mcl[select.cluster]

## combine the sequence objects for all Exomes
EF_120_exons <- falExome[nchar(falExome)>120]

## how many bases are it all toghether again
sum(nchar(EF_120_exons))

## lets select the close to 1:1:1 orthologs
allPreSelectedExons <- EF_120_exons[unlist(Stmcl)]
allPreSelectedExons <- allPreSelectedExons[!is.na(allPreSelectedExons)]

sum(nchar(allPreSelectedExons))

## in contrast to earlier versions I do the bait selection and complete evaluation only on the baits themselves.

SelectedExons <- allPreSelectedExons

############################## START SELECTING BAITS NOT EXONS ##################


contig <- gsub("(EfaB_\\d+):.*", "\\1", names(SelectedExons))
start <- as.numeric(gsub("EfaB_\\d+:(\\d+)-(\\d+)",
                          "\\1", names(SelectedExons)))
end <- as.numeric(gsub("EfaB_\\d+:(\\d+)-(\\d+)",
                       "\\2", names(SelectedExons)))

## just take the central region in each exon and make it a multiple of 120bases long
remainder <- (end-start)%%120
core.start <- start+floor(remainder/2)
core.end <- end-ceiling(remainder/2)

regions.frame <- as.data.frame(cbind(contig, core.start, core.end))

bait.pos <- apply(regions.frame, 1, function (x) {
    along <- seq(as.numeric(x["core.start"]),
                 as.numeric(x["core.end"]), by = 120)
    starts <- along[1:(length(along)-1)]
    ends <- along[2:length(along)]-1
    cbind(contig=x["contig"], starts, ends)
})

bait.pos <- as.data.frame(do.call(rbind, bait.pos))
bait.pos$contig <- as.character(bait.pos$contig)
bait.pos$starts <- as.numeric(as.character(bait.pos$starts))
bait.pos$ends <- as.numeric(as.character(bait.pos$ends))

baits <- unlist(lapply(1:nrow(bait.pos),
                       function(i){
                           substr(falGenome[as.character(bait.pos[i, "contig"])],
                                  as.numeric(bait.pos[i, "starts"]),
                                  as.numeric(bait.pos[i, "ends"]))
                       }))

names(baits) <- paste(bait.pos$contig, ":" ,
                      bait.pos$start, "-",
                      bait.pos$end, sep = "")
table(nchar(baits))

## write initial bait sequences
## write.ESeq(baits, "MYbaits_Eimeria_baits_raw.fasta")

bait.ver <- read.blast6("./baits_vs_ver.blt")
bait.nie <- read.blast6("./baits_vs_nie.blt")
bait.ortho <- rbind(bait.nie, bait.ver)

## ## sort out those with blast hits agains nt:
## ## BLAST screening against all Non-Apicomplexan data in NCBI nt
baits.nt <- read.blast6("baits_vs_nt_wo_Api.blt")
baits <- baits[!names(baits)%in%baits.nt$query_id]
sum(nchar(baits))


## positive selection based on blast agains vermifromis and nieschulzi
baits <- baits[names(baits)%in%bait.ver$query_id]
baits <- baits[names(baits)%in%bait.nie$query_id]
sum(nchar(baits))

### Sort out because of short blast alignments
missed <- by(bait.ortho, bait.ortho$query_id, function (x){
    short <- 120 - min (x[, "alignment_length" ])
})

missed <- abs(unlist(missed))
table(missed>10)
too.much.missed <- names(missed[missed>10])

baits <- baits[!names(baits)%in%too.much.missed]
sum(nchar(baits))


## Sort out based on low percent idendity
lowPer <- by(bait.ortho, bait.ortho$query_id, function (x){
    min (x[, "identity" ])
})

lowPer <- unlist(lowPer)
table(lowPer<88)
too.lowPer <- names(lowPer[lowPer<88])

baits <- baits[!names(baits)%in%too.lowPer]
sum(nchar(baits))

## sort out baits with repeats CAG of GTC
## ## they are already revcom, so only CAG is found
baits <- baits[!grepl("(CAG){3,}|(GTC){3,}" ,baits)]
sum(nchar(baits))

## sort out baits matching the E.faliformis more then once
bait.fal <- read.blast6("./baits_vs_falciformis.blt")
dupes <- bait.fal$query_id[duplicated(bait.fal$query_id)]
baits <- baits[!names(baits)%in%dupes]
sum(nchar(baits))

## select based on gc bigger than 42 but smaller than 65
baits <- baits[ !(get.gc(baits)<42) & !((get.gc(baits)>65)) ]


contig.final <- gsub("(EfaB_\\d+):.*", "\\1", names(baits))
start.final <- as.numeric(gsub("EfaB_\\d+:(\\d+)-(\\d+)",
                          "\\1", names(baits)))
end.final <- as.numeric(gsub("EfaB_\\d+:(\\d+)-(\\d+)",
                       "\\2", names(baits)))


gr.final <- GRanges(seqnames = contig.final, 
                    ranges = IRanges(start=start.final, end = end.final,
                        names = names(baits)),
                    strand = rep("+", times = length(baits)))

## export.gff(gr.final, "MYbaits_Eimeria_V1.single120.gff")

## write.ESeq(baits, "MYbaits_Eimeria_V1.single120.fasta")

final.starts <- gr.final@ranges@start

next.same.ex <- function (x, y) {x + 120 == y}
num.in.ex <- vector()
row <- 1

for (i in 1:(length(final.starts)-1)){
  if(next.same.ex(final.starts[i], final.starts[i+1])){
    row <- row+1
  }
  else {
    num.in.ex[i] <- row
    row <- 1
  }
}

num.in.ex <- as.data.frame(num.in.ex[!is.na(num.in.ex)])
names(num.in.ex) <- "baits"
nrow(num.in.ex)



ggplot(num.in.ex, aes(baits))+
  geom_histogram(binwidth=1, color="white")


gr.final.df <- as.data.frame(gr.final)

marker.dist <- by(gr.final.df, gr.final.df$seqnames, function(x){
  x$start[-1] - x$end[-(nrow(x))]
})

marker.dist <- as.data.frame(cbind(contig=rep(names(marker.dist),
                                        times=unlist(lapply(marker.dist, length))),
                                      between.baits=unlist(marker.dist)))

marker.dist$between.baits <- as.numeric(as.character(marker.dist$between.baits))

top.contigs <- names(tail(table(marker.dist$contig)
                         [order(table(marker.dist$contig))], n=15))

## to have distance only between baits in different exons

ggplot(marker.dist[marker.dist$contig%in%top.contigs,], aes(x=contig, y=between.baits)) +
  geom_boxplot()


