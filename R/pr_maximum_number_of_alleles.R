#' @title Compute the probability distribution of the maximum number of distinct alleles in a DNA mixture
#'
#' @param contributors Character vector with unique names of contributors. Valid names are "U1", "U2", ... for unrelated contributors or the names of pedigree members for related contributors.
#' @param freqs Allele frequencies (see [read_allele_freqs])
#' @param pedigree (optionally) [ped][pedtools::ped] object
#' @param dropout_prs Numeric vector. Dropout probabilities per contributor. Defaults to zeroes.
#' @param fst Numeric. Defaults to 0.
#' @param loci Character vector of locus names (defaults to names attr. of `freqs`)
#' @details A DNA mixture of \eqn{n} contributors contains \eqn{2n} *independent* alleles per locus if the contributors are unrelated; fewer if they are related. This function computes the probability distribution of the maximum number of *distinct* alleles observed across all loci. This is also known as MAC (Maximum Allele Count). Mixture contributors may be related according to an optionally specified pedigree. Optionally, a sub-population correction may be applied by setting `fst>0`.
#'
#' The case where all contributors are unrelated was discussed by Tvedebrink (2014) and is implemented in the `DNAtools` package. Kruijver & Curran (2022)
#' extended this to include related contributors by exploiting the [multiPersonIBD][ribd::multiPersonIBD] function in the `ribd` package.
#'
#' @returns A named numeric vector describing the probability distribution. Numeric values are the probabilities corresponding to the names describing integer values.
#' @examples
#' # define a pedigree of siblings S1 and S2 (and their parents)
#' ped_sibs <- pedtools::nuclearPed(children = c("S1", "S2"))
#'
#' # define allele frequencies
#' freqs <- list(locus1 = c(0.1, 0.9),
#'               locus2 = c(0.25, 0.5, 0.25))
#'
#' # compute dist. of number of alleles for two siblings and one unrelated persons
#' pr_max_number_of_distinct_alleles(contributors = c("S1","S2","U1"), freqs,
#'                                     pedigree = ped_sibs)
#'
#' ## GlobalFiler example (2 unrelated contributors)
#' freqs <- read_allele_freqs(system.file("extdata","FBI_extended_Cauc.csv",
#' package = "numberofalleles"))
#'
#' gf_loci <- c("D3S1358", "vWA", "D16S539", "CSF1PO", "TPOX", "D8S1179",
#'              "D21S11",  "D18S51", "D2S441", "D19S433", "TH01", "FGA",
#'              "D22S1045", "D5S818", "D13S317", "D7S820", "SE33",
#'              "D10S1248", "D1S1656", "D12S391", "D2S1338")
#'
#' p_gf <- pr_max_number_of_distinct_alleles(contributors = c("U1", "U2"),
#'                                             freqs = freqs, loci = gf_loci)
#'
#' barplot(p_gf)
#'
#' @references
#' M. Kruijver & J.Curran (2022). 'The number of alleles in DNA mixtures with related
#' contributors', manuscript submitted
#'
#' T. Tvedebrink (2014). 'On the exact distribution of the number of
#' alleles in DNA mixtures', International Journal of Legal Medicine; 128(3):427--37.
#' \doi{10.1007/s00414-013-0951-3}
#' @export
pr_max_number_of_distinct_alleles <-
  function(contributors, freqs, pedigree,
           dropout_prs = rep(0., length(contributors)),
           fst = 0, loci = names(freqs)){

    pr_distinct_by_locus <- pr_total_number_of_distinct_alleles(contributors = contributors,
                                                                freqs = freqs,
                                                                pedigree = pedigree,
                                                                dropout_prs = dropout_prs,
                                                                fst = fst, loci = loci)$by_locus
    # standardise
    max_number_of_independent_alleles <- 2 * length(contributors)

    p_by_locus <- list()
    cum_p_by_locus <- list()

    for (locus in names(pr_distinct_by_locus)){
      p_locus <- pr_distinct_by_locus[[locus]]

      p_locus_standardised <- setNames(rep(0., max_number_of_independent_alleles + 1),
                                       0:max_number_of_independent_alleles)
      p_locus_standardised[1 + as.integer(names(p_locus))] <- p_locus

      p_by_locus[[locus]] <- p_locus_standardised
      cum_p_by_locus[[locus]] <- cumsum(p_locus_standardised)
    }

    cum_p_across_loci <- Reduce("*", cum_p_by_locus)

    p_mac <- setNames(c(cum_p_across_loci[1], diff(cum_p_across_loci)),
             names(cum_p_across_loci))

    p_mac[p_mac > 0]
}
