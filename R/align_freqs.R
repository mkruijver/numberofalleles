align_freqs_across_populations <- function(freqs_by_population, loci){


  # ensure loci are available in all populations
  for (i_pop in seq_along(freqs_by_population)){
    loci_available <- loci %in% names(freqs_by_population[[i_pop]])

    if (!all(loci_available)){
      stop(paste(names(freqs_by_population[[i_pop]])[!loci_available],
                 collapse = ", "), " not available in population ", i_pop)
    }
  }

  number_of_populations <- length(freqs_by_population)

  aligned_freqs <- list()

  for (locus in loci){
    # collect labels across all pops
    unsorted_allele_labels <- unique(c(lapply(freqs_by_population, function(freqs) names(freqs[[locus]])), recursive = TRUE))
    allele_labels <- sort_allele_labels(unsorted_allele_labels)

    # line them all up in a single matrix
    f <- matrix(data = 0., nrow = length(allele_labels), ncol = number_of_populations,
                dimnames = list(allele_labels, names(freqs_by_population)))

    for (i_pop in seq_len(number_of_populations)){
      f_pop <- freqs_by_population[[i_pop]][[locus]]

      f[match(names(f_pop), allele_labels), i_pop] <- f_pop
    }

    aligned_freqs[[locus]] <- f
  }

  aligned_freqs
}

