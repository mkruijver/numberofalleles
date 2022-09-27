
sample_mixture_alleles_locus_brute_force <- function(freqs_locus, number_of_contributors, number_of_samples){
  g <- numberofalleles:::enumerate_genotypes(number_of_alleles = length(freqs_locus))
  g_pr <- numberofalleles:::get_genotype_probabilities(g = g, f = freqs_locus)

  genotype_samples <- replicate(n = number_of_contributors,
                            g[sample.int(n = length(g_pr), size = number_of_samples, replace = TRUE, prob = g_pr),,drop=FALSE], simplify = FALSE)

  genotype_samples_matrix <- do.call(cbind, genotype_samples)

  allele_labels <- names(freqs_locus)
  allele_levels <- seq_along(allele_labels)

  apply(genotype_samples_matrix, 1, function(x) factor(sort(unique(x))+1L, levels = allele_levels,
                                                       labels = allele_labels), simplify = FALSE)
}

sample_total_allele_count <- function(freqs, number_of_contributors, number_of_samples){

  samples <- lapply(freqs, function(f) sample_allele_count_locus(f, number_of_contributors, number_of_samples))

  Reduce("+", samples)
}

sample_maximum_allele_count <- function(freqs, number_of_contributors, number_of_samples){

  samples <- lapply(freqs, function(f) sample_allele_count_locus(f, number_of_contributors, number_of_samples))

  Reduce(pmax, samples)
}
