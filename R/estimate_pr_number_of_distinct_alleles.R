#' @title Estimate the probability distribution of the number of distinct alleles observed in a DNA mixture for a locus
#'
#' @param number_of_independent_alleles Integer. Number of independent alleles in the mixture.
#' @param f Numeric vector with allele frequencies
#' @param fst Numeric value for sub-population correction (also known as theta)
#' @description For a given number of *independent* alleles, estimates the probability distribution of the number of *distinct* alleles observed in a DNA mixture.
#' @details Due to allele sharing between DNA mixture contributors, the number of *distinct* alleles observed in a mixture is often less than the number of independent alleles in the mixture. For example, if mixture comprises two unrelated contributors, there are four independent alleles. Some of these four independent alleles may be of the same allelic type so that at least one and at most four distinct alleles are observed.
#'
#' This function esimates the probability distribution of the number of *distinct* alleles observed when the mixtures comprises a given number of *independent* alleles. Optionally, a sub-population correction may be applied by setting `fst>0`.
#'
#' The function uses a Monte Carlo estimation procedure.
#' @returns A named numeric vector describing the probability distribution. Numeric values are the probabilities corresponding to the names describing integer values.
#' @export
estimate_pr_number_of_distinct_alleles <- function(number_of_independent_alleles, f,
                                          fst = 0, number_of_samples = 1e5, importance_sampling = FALSE){

  if (!is.numeric(f)){
    stop("f needs to be a numeric vector of allele frequencies")
  }

  if (abs(sum(f) - 1) > 1e-6){
    warning("f does not sum to 1")
  }

  if (length(number_of_independent_alleles) != 1){
    stop("number_of_independent_alleles needs to have length 1")
  }

  if (!(is.numeric(number_of_independent_alleles) | is.integer(number_of_independent_alleles))){
    stop("number_of_independent_alleles needs to be integer valued")
  }

  if (as.character(number_of_independent_alleles) != as.character(as.integer(number_of_independent_alleles))){
    stop("number_of_independent_alleles needs to be integer valued")
  }

  if (number_of_independent_alleles < 0){
    stop("number_of_independent_alleles needs to be non-negative")
  }

  if (length(fst) != 1){
    stop("fst needs to be length 1")
  }

  if (!is.numeric(fst)){
    stop("fst needs to be numeric")
  }

  if (fst > 1){
    stop("fst > 1")
  }

  if (fst < 0){
    stop("fst < 0")
  }

  if (fst != 0){
    stop("fst != 0 is not implemented yet")
  }

  if (number_of_independent_alleles == 0){
    return(stats::setNames(1., "0"))
  }

  if (!importance_sampling){

    return(table(sample_allele_count_num_indep_alleles_locus(freqs_locus = f,
                                                number_of_indep_alleles = number_of_independent_alleles,
                                                number_of_samples = number_of_samples)) / number_of_samples)
  }else{

    pr <- sapply(1:number_of_independent_alleles, function(k){

      m <- numberofalleles:::sample_k_allele_mixtures(f, number_of_alleles = k, number_of_samples = number_of_samples)

      m_pr_n <- pr_mixtures_given_n_independent_alleles(number_of_independent_alleles = number_of_independent_alleles,
                                                          theta = 0, mixtures = m, f = f)

      m_pr_k <- numberofalleles:::pr_k_allele_mixtures_given_n_independent_alleles(number_of_independent_alleles = k, x = m, f = f)

      mean(m_pr_n / m_pr_k)
    })

    return(setNames(pr, 1:number_of_independent_alleles))
  }
}
