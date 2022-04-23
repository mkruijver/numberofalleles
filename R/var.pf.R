#' Variance method for objects of class pf
#'
#' @param x Object of class \code{pf}: output of [pr_total_number_of_distinct_alleles]
#' @param by_locus If \code{TRUE} then the results will be returned locus by locus
#' @param ... other arguments that may p
#'
#' @return Either a vector of variances, one for each locus, or the sum of the locus variances.
#' The variances are the variances of the number of alleles observed at each locus.
#' @rdname var
#'
#' @import methods
#' @seealso [mean.pf]
#' @export
#'
#' @examples
#' freqs <- read_allele_freqs(system.file("extdata","FBI_extended_Cauc.csv",
#'                            package = "numberofalleles"))
#' p <- pr_total_number_of_distinct_alleles(contributors = c("U1","U2"),
#'                                          freqs = freqs)
#' mean(p)
#' var(p)
#'
var <- function(x, ...) {
  if (inherits(x, "pf")) {
    UseMethod('var', x)
  }
  else {
    stats::var(x)
  }
}
#' Variance method for \code{pf} object
#'
#' @rdname var
#' @export
var.pf = function(x, by_locus = FALSE, ...){

  var_by_locus <- sapply(x$by_locus, function(fx) {

    x <- as.numeric(names(fx))
    mean <- sum(fx * x)

    sum(fx * (x - mean)^2)
  })

  if(by_locus){
    return(var_by_locus)
  }else{
    return(sum(var_by_locus))
  }
}
