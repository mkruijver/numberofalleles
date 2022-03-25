#' @title Compute the expectation of the number of distinct alleles observed in a DNA mixture for a locus
#'
#' @param number_of_independent_alleles Integer. Number of independent alleles in the mixture.
#' @param f Numeric vector with allele frequencies
#' @param fst Numeric value for sub-population correction (also known as theta)
#' @description For a given number of *independent* alleles, compute the mean and variance of the number of *distinct* alleles observed in a DNA mixture.
#' @details Due to allele sharing between DNA mixture contributors, the number of *distinct* alleles observed in a mixture is often less than the number of independent alleles in the mixture. For example, if mixture comprises two unrelated contributors, there are four independent alleles. Some of these four independent alleles may be of the same allelic type so that at least one and at most four distinct alleles are observed.
#'
#' This function computes the probability distribution of the number of *distinct* alleles observed when the mixtures comprises a given number of *independent* alleles. Optionally, a sub-population correction may be applied by setting `fst>0`.
#'
#' @examples
#' f <- c(A = 0.1, B = 0.2, C = 0.7)
#'
#' p <- pr_number_of_distinct_alleles(3, f)
#' p_by_hand <- c(sum(f^3), 1 - sum(f^3) - 6 * prod(f), 6 * prod(f))
#' stopifnot(all.equal(as.vector(p), p_by_hand))
#' @export
expected_number_of_distinct_alleles <- function(number_of_independent_alleles,
                                                f, fst = 0){

  moments_check_arguments(number_of_independent_alleles, f, fst)

  if (number_of_independent_alleles == 0){
    return(0)
  }

  p0 <- compute_p0(f, number_of_independent_alleles, fst)
  p1 <- 1 - p0

  sum(p1)
}

NULL
#' @title Compute the variance of the number of distinct alleles observed in a DNA mixture for a locus
#'
#' @param number_of_independent_alleles Integer. Number of independent alleles in the mixture.
#' @param f Numeric vector with allele frequencies
#' @param fst Numeric value for sub-population correction (also known as theta)
#' @description For a given number of *independent* alleles, compute the mean and variance of the number of *distinct* alleles observed in a DNA mixture.
#' @details Due to allele sharing between DNA mixture contributors, the number of *distinct* alleles observed in a mixture is often less than the number of independent alleles in the mixture. For example, if mixture comprises two unrelated contributors, there are four independent alleles. Some of these four independent alleles may be of the same allelic type so that at least one and at most four distinct alleles are observed.
#'
#' This function computes the probability distribution of the number of *distinct* alleles observed when the mixtures comprises a given number of *independent* alleles. Optionally, a sub-population correction may be applied by setting `fst>0`.
#'
#' @examples
#' f <- c(A = 0.1, B = 0.2, C = 0.7)
#' v <- variance_number_of_distinct_alleles(3, f)
#'
#' # now compute variance by hand from the full pr. dist of N
#' p <- pr_number_of_distinct_alleles(3, f)
#'
#' n <- as.numeric(names(p))
#' p_n <- as.vector(p)
#'
#' # compute expected value
#' ev <- expected_number_of_distinct_alleles(3, f)
#'
#' # and variance by hand
#' v2 <- sum(p_n * (n - ev)^2)
#'
#' stopifnot(all.equal(v,v2))
#' @export
variance_number_of_distinct_alleles <- function(number_of_independent_alleles,
                                                f, fst = 0){

  moments_check_arguments(number_of_independent_alleles, f, fst)

  if (number_of_independent_alleles == 0 || fst == 1){
    return(0)
  }

  # determine all indices a<b
  ab <- expand.grid(seq(f), seq(f))
  ab <- ab[ab[,1] < ab[,2],]

  p0 <- compute_p0(f, number_of_independent_alleles, fst)
  p1 <- 1 - p0

  # compute the covariance for each pair
  ab_cov <- apply(ab, 1, function(ab) {
    compute_cov(ab[1], ab[2],
                number_of_independent_alleles,
                p0, p1, fst, f)},
    simplify = TRUE)

  a_var <- p0 * p1

  sum(a_var) + 2 * sum(ab_cov)
}

compute_cov <- function(a, b, m, p0, p1, fst, f){
  compute_p_ia1_ib1(a, b, m, fst, f, p0, p1) - p1[a] * p1[b]
}

compute_p_ia1_ib1 <- function(a, b, m, fst, f, p0, p1){
  j <- 1:m

  1 - (1-compute_p0_a_given_b0(a,b,m,fst,f)) * p0[b] -
      (1-compute_p0_a_given_b0(b,a,m,fst,f)) * p0[a] -
      compute_p0_a(f[a] + f[b], fst, j)
}

compute_p0_a_given_b0 <- function(a,b,m,fst,f){

  if (fst==0){
    return(((1-f[a]-f[b])/(1-f[b]))^m)
  }

  alpha <- (1-fst)/fst
  prod((alpha * (1-f[a]-f[b])+ (0:(m-1) ))/( alpha*(1-f[b]) + (0:(m-1))))
}

compute_p0 <- function(f, m, fst){

  if (fst==1) return(1-f)

  j <- 1:m

  sapply(f, function(f_a){
    compute_p0_a(f_a, fst, j)
  })
}

compute_p0_a <- function(f_a, fst, j){
  prod(((j-1) * fst + (1-fst) * (1-f_a)) / (1 + (j-2) * fst))
}

moments_check_arguments <- function(number_of_independent_alleles,
                            f, fst){
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
}
