#' Mean method for objects of class pf
#'
#' @param x Object of class \code{pf}: output of [pr_total_number_of_distinct_alleles]
#' @param by_locus If \code{TRUE} then the results will be returned locus by locus
#' @param ... other arguments that may p
#'
#' @return either a vector of means, one for each locus, or the sum of the locus means.
#' The means are the expected number of alleles observed at each locus
#' @export
#'
#' @examples
#' freqs <- read_allele_freqs(system.file("extdata","FBI_extended_Cauc.csv",
#'                            package = "numberofalleles"))
#' p <- pr_total_number_of_distinct_alleles(contributors = c("U1","U2"),
#'                                          freqs = freqs)
#' mean(p)
#' var(p)
mean.pf = function(x, by_locus = FALSE, ...){
  mean_by_locus = sapply(x$by_locus, function(loc)sum(loc * as.numeric(names(loc))))

  if(by_locus) mean_by_locus else sum(mean_by_locus)
}
