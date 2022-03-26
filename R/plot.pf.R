#' Plot method for objects of class pf
#'
#' @param x An object of class pf
#' @param ... Any other arguments that should be sent to \codeP
#' @export
plot.pf = function(x, plot.lines = TRUE, plot.points = TRUE, ...){
  xVals = x$noa
  yVals = x$pf
  plot(xVals, yVals, pch = 18, type = "n", ...)
  if(plot.points)
    plot(xVals, yVals, ...)
  if(plot.lines)
  lines(rep(xVals, 2), rep(0, length(yVals)), ...)
}
