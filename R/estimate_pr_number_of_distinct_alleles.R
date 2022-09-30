#' @title Estimate the probability distribution of the number of distinct alleles observed in a DNA mixture for a locus
#'
#' @param number_of_independent_alleles Integer. Number of independent alleles in the mixture.
#' @param f Numeric vector with allele frequencies
#' @param fst Numeric value for sub-population correction (also known as theta)
#' @description For a given number of *independent* alleles, estimates the probability distribution of the number of *distinct* alleles observed in a DNA mixture.
#' @details Due to allele sharing between DNA mixture contributors, the number of *distinct* alleles observed in a mixture is often less than the number of independent alleles in the mixture. For example, if mixture comprises two unrelated contributors, there are four independent alleles. Some of these four independent alleles may be of the same allelic type so that at least one and at most four distinct alleles are observed.
#'
#' This function estimates the probability distribution of the number of *distinct* alleles observed when the mixtures comprises a given number of *independent* alleles. Optionally, a sub-population correction may be applied by setting `fst>0`.
#'
#' The function uses a Monte Carlo estimation procedure.
#' @returns A named numeric vector describing the probability distribution. Numeric values are the probabilities corresponding to the names describing integer values.
#' @examples
#'
#'freqs <- read_allele_freqs(system.file("extdata",
#'                                       "FBI_extended_Cauc.csv",
#'                                       package = "numberofalleles"))
#'
#'# compute the exact distribution of distinct alleles for a 5-person mixture at SE33
#'p <- pr_number_of_distinct_alleles(number_of_independent_alleles = 10,
#'                                   f = freqs$SE33)
#'
#'plot(seq(p), p, type = "b", log = "y", xlab = "Number of distinct alleles")
#'grid()
#'
#'# next we compare it to a simple Monte Carlo estimate and the Importance Sampling estimate
#'set.seed(1)
#'p_hat <- estimate_pr_number_of_distinct_alleles(number_of_independent_alleles = 10,
#'                                                f = freqs$SE33, number_of_samples = 1e4)
#'p_hat_is <- estimate_pr_number_of_distinct_alleles(number_of_independent_alleles = 10,
#'                                                   f = freqs$SE33, number_of_samples = 1e4,
#'                                                   importance_sampling = TRUE)
#'plot(p_hat_is/p-1)
#'# Monte Carlo fails to estimate the tiny probabilities
#'# but Importance Sampling is fairly accurate even with only 10k samples
#'plot(seq(p), p_hat/p, type="b", lty=2, xlab="Number of distinct alleles")
#'lines(p_hat_is/p, type="b")
#'grid()
#'legend("bottomright", legend=c("Monte Carlo", "Importance Sampling"), lty = c(2,1))
#'
#' @export
estimate_pr_number_of_distinct_alleles <- function(number_of_independent_alleles, f,
                                          fst = 0, number_of_samples = 1e5, importance_sampling = FALSE,
                                          uniform_sampling = FALSE){

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

    estimate <- tabulate(sample_allele_count_num_indep_alleles_locus(freqs_locus = f,
                                                                      number_of_indep_alleles = number_of_independent_alleles,
                                                                      number_of_samples = number_of_samples),
                         nbins = number_of_independent_alleles) / number_of_samples

    return(setNames(estimate, seq_len(number_of_independent_alleles)))
  }else{

    pr <- sapply(1:number_of_independent_alleles, function(k){

      number_of_possible_alleles <- sum(f > 0)
      if (number_of_possible_alleles >= k){
        # sample samples of exactly k alleles (without replacement)

        if (uniform_sampling){
          m <- numberofalleles:::sample_k_allele_mixtures_uniform(f, number_of_alleles = k, number_of_samples = number_of_samples)

          m_pr_k <- 1 / choose(length(f), k)
        }else{
          m <- numberofalleles:::sample_k_allele_mixtures(f, number_of_alleles = k, number_of_samples = number_of_samples)

          # and the actual pr. under the sampling scheme
          m_pr_k <- numberofalleles:::pr_k_allele_mixtures(x = m, f = f) * factorial(k)
        }

        # compute the pr. of each sample given n independent alleles
        m_pr_n <- numberofalleles:::pr_mixtures_given_n_independent_alleles(number_of_independent_alleles = number_of_independent_alleles,
                                                            theta = 0, mixtures = m, f = f)

        return(mean(m_pr_n / m_pr_k))
      }
      else if (number_of_possible_alleles < k){
        return(0)
      }

    })

    return(setNames(pr, seq_len(number_of_independent_alleles)))
  }
}
