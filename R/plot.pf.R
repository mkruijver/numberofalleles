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
#' @examples
#'
#' ## Load some data
#' data("FBIfreqs")
#'
#' ## Compute the pedigrees for two and three siblings respectively
#' pedSibs2 = pedtools::nuclearPed(father = "F", mother = "M", children = c("S1", "S2"))
#' pedSibs3 = pedtools::addChildren(father = "F", mother = "M", pedSibs2, ids = "S3")
#'
#' ## Compute the probability functions for each of 2 Unknowns, and 2 Sibs
#' p2U = pr_total_number_of_distinct_alleles(contributors = c("U1", "U2"), freqs = FBIfreqs)
#' p2S =  pr_total_number_of_distinct_alleles(contributors = c("S1", "S2"), freqs = FBIfreqs, pedigree = pedSibs2)
#'
#' ## And plot the two probability distribution functions on top of each other
#' plot(p2U,
#'      line.col = "red",
#'      point.col = "red",
#'      pch = 18,
#'      lwd = 2,
#'      ylim = c(0, 0.15),
#'      xlab = "Total Number of Alleles (TAC)",
#'      ylab = expression(Pr(N == n~"|"~M == m, P)))
#'
#' plot(p2S,
#'      add = TRUE,
#'      point.col = "blue",
#'      line.col = "blue",
#'      lwd = 2)
#'
#' legend("topleft", legend = c("2 U", "2 S"), fill = c("red", "blue"), bty = "n")
#'
#' ## Compute the LR for the number of peaks given the number of contributors and the pedigrees
#' lr = p2U$pf / p2S$pf
#' data.df = data.frame(log10lr = log10(lr), noa = p2U$noa)
#'
#' ## Plot the LR and a grid so that it's easy to see where the LR becomes greater than 1
#' plot(log10lr~noa,
#'      data = data.df,
#'      axes = FALSE,
#'      xlab = "Total number of alleles, n",
#'      ylab = expression(log[10](LR(P[1],P[2]~"|"~N==n, M==m))),
#'      xaxs = "i", yaxs = "i",
#'      xlim = c(22,93),
#'      pch = 18)
#' axis(1)
#' y.exponents = seq(-20, 15, by = 5)
#' for(i in 1:length(y.exponents)){
#'   if(y.exponents[i] == 0)
#'     axis(2, at=y.exponents[i], labels="1", las = 1)
#'   else if(y.exponents[i] == 1)
#'     axis(2, at=y.exponents[i], labels="10", las = 1)
#'   else
#'     axis(2, at = y.exponents[i],
#'          labels = eval(substitute(
#'            expression(10^y),
#'            list(y = y.exponents[i]))),
#'          las = 1)
#' }
#' grid()
#' box()
#'
#' ## Let's look at 2 sibs versus 3 sibs
#'
#' p3S =  pr_total_number_of_distinct_alleles(contributors = c("S1", "S2", "S3"), freqs = FBIfreqs, pedigree = pedSibs3)
#' plot(p3S,
#'      line.col = "green",
#'      point.col = "green",
#'      pch = 18,
#'      lwd = 2,
#'      ylim = c(0, 0.15),
#'      xlab = "Total Number of Alleles (TAC)",
#'      ylab = expression(Pr(N == n~"|"~M == m, P)))
#' plot(p2S,
#'      add = TRUE,
#'      pch = 18,
#'      point.col = "blue",
#'      line.col = "blue",
#'      lwd = 2)
#' legend("topleft", legend = c("3 S", "2 S"), fill = c("green", "blue"), bty = "n")
#'
#' ## And finally two sibs and one unknown versus three sibs
#' p2SU =  pr_total_number_of_distinct_alleles(contributors = c("S1", "S2", "U1"), freqs = FBIfreqs, pedigree = pedSibs2)
#' plot(p3S,
#'      line.col = "green",
#'      point.col = "green",
#'      pch = 18,
#'      lwd = 2,
#'      ylim = c(0, 0.15),
#'      xlab = "Total Number of Alleles (TAC)",
#'      ylab = expression(Pr(N == n~"|"~M == m, P)))
#' plot(p2SU,
#'      add = TRUE,
#'      pch = 18,
#'      point.col = "orange",
#'      line.col = "orange",
#'      lwd = 2)
#' legend("topleft", legend = c("3 S", "2 S + U"), fill = c("green", "orange"), bty = "n")
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
