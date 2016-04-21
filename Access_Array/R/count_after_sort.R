library(pheatmap)
library(reshape)
library(RSvgDevice)


## count before sort first:
## files <- list.files("raw_fastq", pattern = ".*_R1_.*\\.fastq\\.gz", full.names = TRUE)

## counts.list <- list()
## for(f in files) {
##     command <- paste('zcat', f,
##                      '|  grep -c ^@M01')
##     print(command)
##     con <- pipe (command)
##     n <- sub("^raw_fastq/(.*-chip\\d).*", "\\1", f) 
##     counts.list[[n]] <- readLines(con)
##     close(con)
## }

## foo <- melt(counts.list)

## write.csv(foo, "/home/ele/Dropbox/Eimeria_Wild/Access_array/raw_count_table.csv",
##           row.names=FALSE)

## Counts after primer matching
files <- list.files("qual_filter/sorted_amps_trimmed/", "*_R1_*.fastq", full.names = TRUE)

counts.list <- list()

for(f in files) {
    command <- paste('grep \'@RECORD\'', f,
                     '| awk -F\':|_R1\' \'{print $2}\' | sort | uniq -c')
    print(command)
    con <- pipe (command)
    counts.list[[f]] <- readLines(con)
    close(con)
}

foo <- melt(counts.list)

foo$value <- sub("^\\s+|\\s+$", "", foo$value)

foo$count <- as.numeric(sapply(strsplit(foo$value, " "), "[", 1))
foo$sample <- as.character(sapply(strsplit(foo$value, " "), "[", 2))

foo$amplicon <- as.character(sapply(foo$L1,  function (x){
    s <- strsplit(x, "/|_")[[1]][c(7, 8)]
    paste(s, collapse = "")
}))

foo$value <- NULL
foo$L1 <- NULL


read.matrix <- reshape(foo,
                       direction = "wide",
                       v.names = "count", timevar = "sample",
                       idvar = "amplicon")

rownames(read.matrix) <- read.matrix$amplicon
read.matrix$amplicon <- NULL
read.matrix[is.na(read.matrix)] <- 0
names(read.matrix) <- gsub("count.", "", names(read.matrix))


oldnames <- names(read.matrix)
chipno <- gsub(".*(chip\\d)", "\\1", names(read.matrix))
names(read.matrix) <- gsub("-chip\\d", "", names(read.matrix))

devSVG("/home/ele/Dropbox/Eimeria_Wild/Access_array/figures_heatmaps/trimmed_chip1_heat.svg")
pheatmap(log10(read.matrix[, chipno%in%"chip1"]+1))
dev.off()

devSVG("/home/ele/Dropbox/Eimeria_Wild/Access_array/figures_heatmaps/trimmed_chip2_heat.svg")
pheatmap(log10(read.matrix[, chipno%in%"chip2"]+1))
dev.off()

devSVG(
       width=20, height = 10)
pheatmap(log10(read.matrix[, chipno%in%c("chip1","chip2")]+1),
         labels_col = oldnames[chipno%in%c("chip1","chip2")])
dev.off()

devSVG("/home/ele/Dropbox/Eimeria_Wild/Access_array/figures_heatmaps/trimmed_chip3_heat.svg")
pheatmap(log10(read.matrix[, chipno%in%"chip3"]+1))
dev.off()

devSVG("/home/ele/Dropbox/Eimeria_Wild/Access_array/figures_heatmaps/trimmed_chip4_heat.svg")
pheatmap(log10(read.matrix[, chipno%in%"chip4"]+1))
dev.off()

devSVG("/home/ele/Dropbox/Eimeria_Wild/Access_array/figures_heatmaps/trimmed_chip3u4_heat.svg",
       width=20, height = 10)
pheatmap(log10(read.matrix[, chipno%in%c("chip3","chip4")]+1),
         labels_col = oldnames[chipno%in%c("chip3","chip4")])
dev.off()


