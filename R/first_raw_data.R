loc <- read.csv("./2015_only_table_final.csv")

A_ids <- as.character(loc$Prot_A_IDs)
B_ids <- as.character(loc$Prot_B_IDs)

A_ids <- strsplit(A_ids, ";")
B_ids <- strsplit(B_ids, ";")

A_ids_N <- unlist(lapply(A_ids, length))
B_ids_N <- unlist(lapply(B_ids, length))

A_ids_C <- unlist(A_ids)
B_ids_C <- unlist(B_ids)

A_reps <- rep(1:nrow(loc), times=A_ids_N)
B_reps <- rep(1:nrow(loc), times=B_ids_N)


mice <- loc[c(A_reps, B_reps), c(1:5, 13, 16) ]
mice <- cbind(c(A_ids_C, B_ids_C), mice)
names(mice)[1] <- "ID_CZ_DB"

## write.csv(mice, "mouse_table_BR_2015.csv", row.names=FALSE)
