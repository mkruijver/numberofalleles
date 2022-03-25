#' @title Compute the probability distribution of the number of independent alleles in a mixture with dropout
#'
#' @param dropout_prs Numeric vector. Dropout probabilities per contributor.
#' @description Without dropout, each mixture contributor has two *independent* but not necessarily *distinct* alleles that are represented in the DNA mixture. If the probability of dropout for a mixture contributor is greater than zero, then the mixture contributor has either 0 (full dropout), 1 (partial dropout) or 2 (no dropout) independent alleles that are represented in the mixture. This function computes the probability distribution of the number of independent alleles that unrelated mixture contributors have in total for a locus given their dropout parameters.
#' @seealso [pr_independent_alleles_ped]
#' @examples
#' # a dropout pr. of 0.5
#' p <- pr_independent_alleles(0.5)
#' stopifnot(all.equal(as.vector(p),
#'                     c(0.5^2, 2 * 0.5 * 0.5, (1-0.5)^2)))
#'
#' # one contrib. without dropout and one with d=0.5
#' p1 <- pr_independent_alleles(c(0, 0.5))
#' stopifnot(identical(as.integer(names(p1)),
#'                     as.integer(names(p)) + 2L))
#' @export
pr_independent_alleles <- function(dropout_prs){

  if (!is.numeric(dropout_prs)){
    stop("dropout_prs needs to be a numeric vector")
  }

  if (any(dropout_prs < 0)){
    stop("dropout_prs needs to be non-negative")
  }

  if (any(dropout_prs > 1)){
    stop("dropout_prs should not exceed 1")
  }

  number_of_non_dropped_alleles_by_u <-
    lapply(dropout_prs, function(d){
      u_with_drop <- prepare_drop_IBD(multi_person_IBD = data.frame(U="1 2"),
                                      person_name = "U", dropout_pr = d)[[1]]

      x <- sapply(u_with_drop$U_retained, length)
      fx <- u_with_drop$U_drop_pr

      u_x <- unique(x)
      u_fx <- sapply(u_x, function(x0) sum(fx[x==x0]))

      stats::setNames(u_fx[u_fx > 0], u_x[u_fx > 0])
    })

  compute_pr_of_sum(number_of_non_dropped_alleles_by_u)
}
