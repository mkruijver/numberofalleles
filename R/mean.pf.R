#' Mean method for objects of class pf
#'
#' @param x Object of class \code{pf}: output of [pr_total_number_of_distinct_alleles]
#' @param bylocus If \code{TRUE} then the results will be returned locus by locus
#' @param ... other arguments that may p
#'
#' @return either a vector of means, one for each locus, or the sum of the locus means.
#' The means are the expected number of alleles observed at each locus
#' @export
#'
#' @examples
#' data(FBIfreqs)
#' p <- pr_total_number_of_distinct_alleles(contributors = c("U1","U2"),
#'                                          freqs = FBIfreqs)
#' mean(p)
#' var(p)
mean.pf = function(x, bylocus = FALSE, ...){
  locusMeans = sapply(x$byLocus, function(loc)sum(loc * as.numeric(names(loc))))

  if(bylocus)
    return(locusMeans)
  ##else
  return(sum(locusMeans))
}
