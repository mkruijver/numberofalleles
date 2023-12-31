% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pr_independent_alleles_ped.R
\name{pr_independent_alleles_ped}
\alias{pr_independent_alleles_ped}
\title{Compute the probability distribution of the number of independent alleles in a mixture with dropout and related contributors}
\usage{
pr_independent_alleles_ped(pedigree, ped_contributors, dropout_prs)
}
\arguments{
\item{pedigree}{\link[pedtools:ped]{ped} object}

\item{ped_contributors}{Character vector with unique names of contributors. Valid names are the names of pedigree members.}

\item{dropout_prs}{Numeric vector. Dropout probabilities per contributor.}
}
\value{
A named numeric vector describing the probability distribution. Numeric values are the probabilities corresponding to the names describing integer values.
}
\description{
When mixture contributors are related according to a pedigree, they may share some alleles identical by descent so that their total number of \emph{independent} alleles is smaller than two times the number of contributors. The number of \emph{independent} alleles can be further reduced if dropout plays a role. This function computes the probability distribution of the number of independent alleles that related mixture contributors have in total for a locus given their dropout parameters. Note that the number of \emph{distinct} alleles that is observed at the locus is typically smaller than the number of independent alleles due to allele sharing.
}
\examples{
# without dropout, a father-mother-child mixture has 4 indep. alleles
p <- pr_independent_alleles_ped(pedtools::nuclearPed(),
                                ped_contributors = as.character(1:3),
                                dropout_prs = rep(0, 3))
stopifnot(identical(p,
                    stats::setNames(1., "4")))
}
\seealso{
\link{pr_independent_alleles}
}
