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
#' data(FBIfreqs)
#' p <- pr_total_number_of_distinct_alleles(contributors = c("U1","U2"),
#'                                          freqs = FBIfreqs)
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

  var_by_locus <- sapply(x$byLocus, function(fx) {

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