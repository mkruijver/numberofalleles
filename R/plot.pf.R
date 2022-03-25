#' Plot method for objects of class pf
#'
#' @param x An object of class pf
#' @param ... Any other arguments that should be sent to barplot
#' @export
plot.pf = function(x, ...){
  barplot(x$pf, names.arg = as.character(x$noa), ...)
  box()
}
