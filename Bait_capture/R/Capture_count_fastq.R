library(reshape)


## count what is in the files
F.filex <- list.files("/SAN/BegenDiv_raw/160304_2015-34-EH-parasites/Data/Intensities/BaseCalls/",
                      "R1_001.fastq.gz", full.names = TRUE)

read.counts.l <- list()

for(file in F.filex){
    command <- paste("zcat", file, "| wc -l")
    read.counts.l[file] <- readLines(pipe(command))
}


read.counts <- lapply(read.counts.l, as.numeric)

read.counts <- lapply(read.counts, "/", 4)

## remove the folder name and the uninformative ending of the filename
names(read.counts) <- 
    gsub("/SAN/BegenDiv_raw/160304_2015-34-EH-parasites/Data/Intensities/BaseCalls//",
         "", names(read.counts))

names(read.counts) <- 
    gsub(".fastq.gz",
         "", names(read.counts))


rc <- melt(read.counts)

## analyse proportions of Eimeria reads reads


## count what is in the files
B.filex <- list.files("/SAN/BegenDiv_raw_first_analysis/sample_subsets/",
                      ".fasta", full.names = TRUE)

## remove the one wired file which does not give any hit this would
## otherwise produce an error in the pipe command below
B.filex <- 
    B.filex[!B.filex%in%"/SAN/BegenDiv_raw_first_analysis/sample_subsets//2812Single_S10_R1_001.fasta"]


read.blast6 <- function (file){
    blast <- read.delim(file, header=FALSE)
    names(blast) <- c("query_id", "subject_id", "identity",
                      "alignment_length", "mismatches",
                      "gap_opens", "q_start", "q_end",
                      "s_start", "s_end", "evalue", "bit_score")
    return(blast)
}

blast.results.l <- list()

for(file in B.filex){
    command <- paste("blastn -query", file,
                     "-db /data/db/blastdb/Eimeria_falciformis/Eimeria_contigs_final.fa" ,
                     "-outfmt 6 -max_target_seqs 1 -num_threads 18")
    print(command)
    blast.results.l[[file]] <- read.blast6(pipe(command))
}

names(blast.results.l) <- 
    gsub("/SAN/BegenDiv_raw_first_analysis/sample_subsets//","",  names(blast.results.l))

names(blast.results.l) <- 
    gsub(".fasta",
         "", names(blast.results.l))


uniq.query <- lapply(blast.results.l, function(x){
                         length(unique(x$query_id))
                     })

uniq.query <- melt(uniq.query)


S.table <- merge(rc,uniq.query, by="L1", all=TRUE)
names(S.table) <- c("library", "total.reads", "match.in.10k")
S.table[is.na(S.table)] <- 0

S.table$expected.overall <- round((S.table$match.in.10k/10000)*S.table$total.reads)

S.table[S.table$total.reads<10000, "expected.overall"] <-
    S.table[S.table$total.reads<10000, "match.in.10k"] 



uniq.subject <- lapply(blast.results.l, function(x){
                           length(unique(x$subject))
                       })

uniq.subject <- melt(uniq.subject)


