#' Plot method for objects of class pf
#'
#' @param x An object of class pf
#' @param plot.lines if \code{TRUE} then lines will be drawn from 0 to the
#' probability value for each x value.
#' @param plot.points if \code{TRUE} then points will be plotted at the
#' probability value for each x value.
#' @param line.col the colour of the lines drawn for each probability mass
#' @param point.col the colour of the points plotted for each probability mass
#' @param ... Any other arguments that should be sent to \code{plot},
#' \code{arrows}, or \code{points}.
#' @param add If \code{TRUE} then the plotting information will be added to
#' the existing plot.
#' @export
plot.pf = function(x,
                   plot.lines = TRUE,
                   plot.points = TRUE,
                   line.col = "black",
                   point.col = "black",
                   add = FALSE,
                   ...){
  xVals = x$noa
  yVals = x$pf

  if(!add)
    plot(xVals, yVals, type = "n", ...)

  if(plot.points)
    points(xVals, yVals,
           col = point.col,
           ...)

  if(plot.lines)
    arrows(x0 = xVals,
           x1 = xVals,
           y0 = rep(0, length(yVals)),
           y1 = yVals,
           length = 0, ## No head
           col = line.col,
           ...
    )
}
