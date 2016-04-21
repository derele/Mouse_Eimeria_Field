library(pheatmap)
library(reshape)
library(RSvgDevice)


## Counts after primer matching
files <- list.files("qual_filter/sorted_amps_trimmed", "*_R1_*.fastq", full.names = TRUE)

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
    s <- strsplit(x, "/|_")[[1]][c(6,7)]
    paste(s, collapse = "")
}))

foo$value <- NULL
foo$L1 <- NULL

foo.sample <- as.data.frame(tapply(as.numeric(as.character(foo$count)),
                                   foo$sample, sum))

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
names(read.matrix) <- gsub("_chip\\d", "", names(read.matrix))

devSVG("/home/ele/Dropbox/Eimeria_Wild/Access_array/figures_heatmaps/chip1_QUALF_heat.svg")
pheatmap(log10(read.matrix[, chipno%in%"chip1"]+1))
dev.off()

devSVG("/home/ele/Dropbox/Eimeria_Wild/Access_array/figures/chip2_QUALF_heat.svg")
pheatmap(log10(read.matrix[, chipno%in%"chip2"]+1))
dev.off()

devSVG("/home/ele/Dropbox/Eimeria_Wild/Access_array/figures/chip1u2_QUALF_heat.svg",
       width=20, height = 10)
pheatmap(log10(read.matrix[, chipno%in%c("chip1","chip2")]+1),
         labels_col = oldnames[chipno%in%c("chip1","chip2")])
dev.off()

devSVG("/home/ele/Dropbox/Eimeria_Wild/Access_array/figures/chip3_QUALF_heat.svg")
pheatmap(log10(read.matrix[, chipno%in%"chip3"]+1))
dev.off()

devSVG("/home/ele/Dropbox/Eimeria_Wild/Access_array/figures/chip4_QUALF_heat.svg")
pheatmap(log10(read.matrix[, chipno%in%"chip4"]+1))
dev.off()

devSVG("/home/ele/Dropbox/Eimeria_Wild/Access_array/figures/chip3u4_QUALF_heat.svg",
       width=20, height = 10)
pheatmap(log10(read.matrix[, chipno%in%c("chip3","chip4")]+1),
         labels_col = oldnames[chipno%in%c("chip3","chip4")])
dev.off()


