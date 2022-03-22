#' @title Compute the probability distribution of the number of distinct alleles observed in a DNA mixture for a locus
#'
#' @param number_of_independent_alleles Integer. Number of independent alleles in the mixture.
#' @param f Numeric vector with allele frequencies
#' @param fst Numeric value for sub-population correction (also known as theta)
#' @param brute_force Logical. Should a brute force algorithm be used?
#' @description For a given number of \emph{independent} alleles, compute the probability distribution of the number of \emph{distinct} alleles observed in a DNA mixture.
#' @details Due to allele sharing between DNA mixture contributors, the number of \emph{distinct} alleles observed in a mixture is often less than the number of independent alleles in the mixture. For example, if mixture comprises two unrelated contributors, there are four independent alleles. Some of these four independent alleles may be of the same allelic type so that at least one and at most four distinct alleles are observed.
#'
#' This function computes the probability distribution of the number of \emph{distinct} alleles observed when the mixtures comprises a given number of \emph{independent} alleles. Optionally, a sub-population correction may be applied by setting \code{fst>0}.
#'
#' An efficient way of computing the probability distribution was given by Tvedebrink (2014) and was slightly adapted by Kruijver & Curran (2022) to handle the case of an odd number of independent alleles. A much slower brute force algorithm is also implemented (argument \code{brute_force=TRUE}) for testing purposes.
#' @examples
#' f <- c(A = 0.1, B = 0.2, C = 0.7)
#'
#' p <- pr_number_of_distinct_alleles(3, f)
#' p_by_hand <- c(sum(f^3), 1 - sum(f^3) - 6 * prod(f), 6 * prod(f))
#' stopifnot(all.equal(as.vector(p), p_by_hand))
#' @references
#' M. Kruijver & J.Curran (2022). 'The number of alleles in DNA mixtures with related
#' contributors', manuscript submitted
#'
#' T. Tvedebrink (2014). 'On the exact distribution of the number of
#' alleles in DNA mixtures', International Journal of Legal Medicine; 128(3):427--37.
#' <https://doi.org/10.1007/s00414-013-0951-3>
#' @export
pr_number_of_distinct_alleles <- function(number_of_independent_alleles, f,
                                          fst = 0,
                                          brute_force = FALSE){

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

  if (number_of_independent_alleles == 0){
    return(stats::setNames(1., "0"))
  }

  partitions <- partitions::parts(number_of_independent_alleles)
  partitions_list <- apply(partitions, 2, function(p) p[p>0], simplify = FALSE)

  weights <- sapply(partitions_list, weights_cpp)

  pr_by_partition <- numeric(length(partitions_list))

  pr_distinct <- stats::setNames(numeric(number_of_independent_alleles),
                                 seq(number_of_independent_alleles))

  for (i in seq_along(partitions_list)){
    alpha <- partitions_list[[i]]

    s <- if (brute_force) S_brute_force(f, alpha, fst) else S_recursive(f, alpha, fst)
    pr_by_partition[i] <- weights[i] * s

    pr_distinct[length(alpha)] <- pr_distinct[length(alpha)] + pr_by_partition[i]
  }

  return(pr_distinct[pr_distinct>0])
}
