seq <- read.ESeq("/data/RAW_SEQ/2015_16_2015_18/qual_filter/sorted_amps_trimmed/flashed/Eim_Nuc9409_F_Eim_Nuc9410_R_flashed.fasta")
seq <- seq[!grepl("foo", names(seq))]

samples <- gsub("RECORD:(.*?_chip\\d)_.*", "\\1", names(seq))

non.uniq <- duplicated(cbind(samples, seq))

samples.dup <- samples[non.uniq]
seq.dup <- seq[non.uniq]

table(samples.dup)

write.ESeq(seq.dup[samples.dup%in%"21898AF_chip4"], "Eim_Nuc9409_sample_21898AF_chip4.fasta")
write.ESeq(seq.dup[samples.dup%in%"37_SK2809_chip2"], "Eim_Nuc9409_sample_SK2808_chip2.fasta")


## for some Ap5 data

seq <- read.ESeq("/data/RAW_SEQ/2015_16_2015_18/qual_filter/sorted_amps_trimmed/flashed/Eim_Ap5_F_Eim_Ap5_R_flashed.fasta")
seq <- seq[!grepl("foo", names(seq))]

seq <- seq[!grepl("foo", names(seq))]

samples <- gsub("RECORD:(.*?_chip\\d)_.*", "\\1", names(seq))

non.uniq <- duplicated(cbind(samples, seq))

samples.dup <- samples[non.uniq]
seq.dup <- seq[non.uniq]

table(samples.dup)

write.ESeq(seq.dup[samples.dup%in%"21898AF_chip4"], "Eim_Ap5_sample_21898AF_chip4.fasta")
write.ESeq(seq.dup[samples.dup%in%"37_SK2809_chip2"], "Eim_Ap5_sample_SK2809_chip2.fasta")

write.ESeq(seq.dup[samples.dup%in%"36_SK2808_chip2"], "Eim_Ap5_sample_SK2808_chip2.fasta")


write.ESeq(seq.dup[samples.dup%in%"16_SK2695_chip1"], "Eim_Ap5_sample_SK2695_chip1.fasta")

### some Mus data

seq <- read.ESeq("/data/RAW_SEQ/2015_16_2015_18/qual_filter/sorted_amps_trimmed/flashed/_Ap5_R_flashed.fasta")
seq <- seq[!grepl("foo", names(seq))]

seq <- seq[!grepl("foo", names(seq))]

samples <- gsub("RECORD:(.*?_chip\\d)_.*", "\\1", names(seq))

non.uniq <- duplicated(cbind(samples, seq))

samples.dup <- samples[non.uniq]
seq.dup <- seq[non.uniq]

table(samples.dup)




