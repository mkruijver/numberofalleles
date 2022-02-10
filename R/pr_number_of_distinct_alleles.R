#' @title Compute probability distribution of number of distinct alleles in a DNA mixture at a locus
#'
#' @param number_of_alleles Integer. Number of (not necessarily distinct) alleles in the mixture.
#' @param f Numeric Vector with allele frequencies
#' @param brute_force Logical Should a brute force algorithm be used?
#' @details TODO
#' @examples
#' f <- c(A = 0.1, B = 0.2, C = 0.7)
#'
#' pr_number_of_distinct_alleles(3, f)
#' @export
pr_number_of_distinct_alleles <- function(number_of_alleles, f,
                                          brute_force = FALSE){

  if (!is.numeric(f)){
    stop("f needs to be a numeric vector of allele frequencies")
  }

  if (abs(sum(f) - 1) > 1e-6){
    warning("f does not sum to 1")
  }

  if (length(number_of_alleles) != 1){
    stop("number_of_alleles needs to have length 1")
  }

  if (!(is.numeric(number_of_alleles) | is.integer(number_of_alleles))){
    stop("number_of_alleles needs to be integer valued")
  }

  if (as.character(number_of_alleles) != as.character(as.integer(number_of_alleles))){
    stop("number_of_alleles needs to be integer valued")
  }

  partitions <- partitions::parts(number_of_alleles)
  partitions_list <- apply(partitions, 2, function(p) p[p>0], simplify = FALSE)

  weights <- sapply(partitions_list, get_weights)

  pr_by_partition <- numeric(length(partitions_list))

  pr_distinct <- stats::setNames(numeric(number_of_alleles),
                                 seq(number_of_alleles))

  for (i in seq_along(partitions_list)){
    alpha <- partitions_list[[i]]

    s <- if (brute_force) S_brute_force(f, alpha) else S_recursive(f, alpha)
    pr_by_partition[i] <- weights[i] * s

    pr_distinct[length(alpha)] <- pr_distinct[length(alpha)] + pr_by_partition[i]
  }

  return(pr_distinct[pr_distinct>0])
}

get_weights <- function(x){
  x <- x[x>0]

  factorial(sum(x)) /
    (prod(factorial(x)) * prod(factorial(rle(sort(x))$lengths)))
}
