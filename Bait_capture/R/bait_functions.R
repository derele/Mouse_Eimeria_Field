## Sorry uncommented code for those who venture here:

read.ESeq <- function (file){
  content <- readLines(file)
  header.lines <- grep( "^>", content)
  start.lines <- header.lines+1
  end.lines <- c(header.lines[-1]-1, length(content))
  sq <- sapply(1:length(start.lines), function (i) {
    list(content[start.lines[i]:end.lines[i]])})
  names(sq) <- substr(content[header.lines],2, nchar(content[header.lines]))
  sq <- unlist(lapply(sq, paste, collapse=""))
  sq <- sq[nchar(sq)>0]
  class(sq) <- "ESeq"
  return(sq)
}

write.ESeq <- function (ESeq.obj,
                        path=paste(getwd(), "/", substitute(ESeq.obj), ".fasta", sep="")){
  write(paste(">",
              names(ESeq.obj), "\n",
              ESeq.obj, sep=""),
        file=path)
}

print.ESeq <- function(x, n=10){
  cat("* ESeq object, easy biological sequence  **\n")
  pfun <- function (x) { if(nchar(x)<10) {
    cat (x, "\n")}
  else{cat(substr(x, 1, 4), "...", "further", nchar(x)-8, "letters ...",
           substr(x, nchar(x)-3, nchar(x)),"\n")}}
  if(length(x)<n){
    sapply (x, pfun)}
  else {
    sapply(x[1:4], pfun)
    cat("...\nfurther", length(x)-8, "sequences\n...\n")
    sapply(x[(length(x)-4):length(x)], pfun)}
}

get.gc <- function (x) (nchar(gsub("A|T|N|a|t|n", "",  x))/
                        nchar(gsub("N|n", "", (x)))*100)

strReverse <- function(y) sapply(lapply(strsplit(y, NULL), rev), paste,
                                 collapse="")

revcom <- function(x){
  r <- strReverse(x)
  return(chartr("acgtryswkmbdhvnxACGTRYSWKMBDHVNX", "tgcayrswmkvhdbnxTGCAYRSWMKVHDBNX", r))
}


read.blast6 <- function (file){
    blast <- read.delim(file, header=FALSE)
    names(blast) <- c("query_id", "subject_id", "identity",
                      "alignment_length", "mismatches",
                      "gap_opens", "q_start", "q_end",
                      "s_start", "s_end", "evalue", "bit_score")
    return(blast)
}

get.Start.End.blast <- function (blast.table, genome) {
    select <- apply(blast.table, 1, function (x) {
        subject <- as.character(x["subject_id"])
        ##  to get the coordinates for the exons in the two
        start <- as.numeric(min(x["s_start"], x["s_end"]))
        end <- as.numeric(max(x["s_start"], x["s_end"]))
        seq <- substring(genome[subject], start, end)
        if(x["s_start"] > x["s_end"]){
          seq <- revcom(seq)
        }
        return(seq)
    })
    names(select) <- paste(blast.table$subject_id, ":",
                           blast.table$s_start, "-",
                           blast.table$s_end, sep = "")
    return(select)
}


get.Start.End.blast.query <- function (blast.table, Exome) {
    select <- apply(blast.table, 1, function (x) {
        query <- as.character(x["query_id"])
        ##  to get the coordinates for the exons in the two
        start <- as.numeric(min(x["q_start"], x["q_end"]))
        end <- as.numeric(max(x["q_start"], x["q_end"]))
        substring(Exome[query], start, end)        
    })
    names(select) <- paste(blast.table$query_id, ":",
                           blast.table$q_start, "-",
                           blast.table$q_end, sep = "")
    return(select)
}


tabulate.mcl <- function (T) {
    T.tab <- lapply(T, function (x) table(substr(x, 1, 4)))
    dat <- data.frame()
    for(i in seq(along=T.tab)) for(j in names(T.tab[[i]])){
        dat[i,j] <- T.tab[[i]][j]}
#    dat[is.na(dat)] <- 0
    return(dat)
}


