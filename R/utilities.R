
get_number_of_unique_alleles <- function(...) length(unique(c(...)))

columns <- function(x) lapply(seq_len(ncol(x)), function(i) x[,i])
