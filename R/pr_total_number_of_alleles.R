#' @title Compute the probability distribution of total number of distinct alleles in a DNA mixture
#'
#' @param contributors Character vector with unique names of contributors. Valid names are "U1", "U2", ... for unrelated contributors or the names of pedigree members for related contributors.
#' @param freqs Allele frequencies (see [read_allele_freqs])
#' @param pedigree (optionally) [ped][pedtools::ped] object
#' @param dropout_prs Numeric vector. Dropout probabilities per contributor. Defaults to zeroes.
#' @param fst Numeric. Defaults to 0.
#' @param loci Character vector of locus names (defaults to names attr. of `freqs`)
#' @details A DNA mixture of \eqn{n} contributors contains \eqn{2n} *independent* alleles per locus if the contributors are unrelated; fewer if they are related. This function computes the probability distribution of the total number of *distinct* alleles observed across all loci. Mixture contributors may be related according to an optionally specified pedigree. Optionally, a sub-population correction may be applied by setting `fst>0`.
#'
#' The case where all contributors are unrelated was discussed by Tvedebrink (2014) and is implemented in the `DNAtools` package. Kruijver & Curran (2022)
#' extended this to include related contributors by exploiting the [multiPersonIBD][ribd::multiPersonIBD] function in the `ribd` package.
#'
#' @examples
#' # define a pedigree of siblings S1 and S2 (and their parents)
#' ped_sibs <- pedtools::nuclearPed(children = c("S1", "S2"))
#'
#' # define allele frequencies
#' freqs <- list(locus1 = c(0.1, 0.9),
#'               locus2 = c(0.25, 0.5, 0.25))
#'
#' # compute dist. of number of alleles for two siblings and one unrelated persons
#' pr_total_number_of_distinct_alleles(contributors = c("S1","S2","U1"), freqs,
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
#' p_gf <- pr_total_number_of_distinct_alleles(contributors = c("U1", "U2"),
#'                                             freqs = freqs, loci = gf_loci)
#'
#' barplot(p_gf)
#' @references
#' M. Kruijver & J.Curran (2022). 'The number of alleles in DNA mixtures with related
#' contributors', manuscript submitted
#'
#' T. Tvedebrink (2014). 'On the exact distribution of the number of
#' alleles in DNA mixtures', International Journal of Legal Medicine; 128(3):427--37.
#' \doi{10.1007/s00414-013-0951-3}
#' @export
pr_total_number_of_distinct_alleles <-
  function(contributors, freqs, pedigree,
           dropout_prs = rep(0., length(contributors)),
           fst = 0, loci = names(freqs)){

  if (!is.character(contributors)){
    stop("contributors should be a character vector")
  }

  if (!missing(pedigree)){
    if (!inherits(pedigree, "ped")){
      stop("pedigree should be of class ped")
    }

    ped_names <- pedigree$ID

    forbidden_names <- paste0("U", seq_len(length(contributors)))

    forbidden_names_in_ped <- intersect(ped_names, forbidden_names)

    if (length(forbidden_names_in_ped) > 0){
      stop("Pedigree contains illegal name(s): ",
           paste(forbidden_names_in_ped, collapse = ", "))
    }

  }else{
    ped_names <- character()
  }

  # verify that names of other contributors are U1, U2, ..
  unr_names <- contributors[!(contributors %in% ped_names)]

  if (length(unr_names) > 0){
    expected_unr_names <- paste0("U", seq_along(unr_names))

    if (!setequal(expected_unr_names, unr_names)){
      stop("Expected unrelated contributor(s) named ",
           paste(expected_unr_names, collapse = ", "),
           " instead of ", paste(unr_names, collapse = ", ") )
    }
  }

  if (!is.list(freqs)){
    stop("freqs should be a list")
  }
  if (!all(sapply(freqs, is.numeric))){
    stop("freqs should be a list of numeric vectors")
  }
  for (locus in loci){
    if (!(locus %in% names(freqs))){
      stop(paste0("freqs not available for locus "), locus)
    }
  }

  if (!is.numeric(dropout_prs)){
    stop("dropout_prs needs to be a numeric vector")
  }

  if (any(dropout_prs < 0)){
    stop("dropout_prs needs to be non-negative")
  }

  if (any(dropout_prs > 1)){
    stop("dropout_prs should not exceed 1")
  }

  if (length(dropout_prs) != length(contributors)){
    stop("contributors and dropout_prs need to have the same length")
  }


  if (!isTRUE(all.equal(rep(1.0, length(loci)),
                        unname(sapply(freqs[loci], sum))))){
    warning("freqs do not sum to 1 at all selected loci")
  }

  # first determine pr. dist. of number of independent alleles
  # for pedigree members and unrelated persons

  ped_contributors <- contributors[contributors %in% ped_names]
  number_of_contributors <- length(contributors)

  if (length(ped_contributors) > 0){

    ped_dropout_prs <- dropout_prs[match(ped_contributors, contributors)]

    pr_num_indep_alleles_ped <- pr_independent_alleles_ped(
        pedigree = pedigree, ped_contributors = ped_contributors,
        dropout_prs = ped_dropout_prs)
  }else{
    # create dummy
    pr_num_indep_alleles_ped <- stats::setNames(1., "0")
  }

  number_of_unrelated_contributors <- length(unr_names)

  if (number_of_unrelated_contributors > 0){
    unr_dropout_prs <- dropout_prs[match(unr_names, contributors)]
    pr_num_indep_alleles_unr <- pr_independent_alleles(unr_dropout_prs)
  }else{
    # create dummy
    pr_num_indep_alleles_unr <- stats::setNames(1., "0")
  }

  # number of indep. alleles is the sum of those for pedigree members
  # and the unrelateds
  pr_num_indep_alleles <- compute_pr_of_sum(
    list(pr_num_indep_alleles_unr, pr_num_indep_alleles_ped))

  pr_by_locus <- list()

  # compute pr. dist. of # distinct alleles locus-wise
  for (locus in loci){
    pr_locus <- stats::setNames(rep(0, 2 * number_of_contributors + 1),
                                0:(2 * number_of_contributors))

    f <- freqs[[locus]]

    for (i in seq_along(pr_num_indep_alleles)){
      num_indep_alleles <- as.integer(names(pr_num_indep_alleles)[i])
      pr <- pr_num_indep_alleles[i]

      pr_distinct <- pr_number_of_distinct_alleles(num_indep_alleles, f, fst = fst)

      idx <- match(names(pr_distinct), names(pr_locus))
      pr_locus[idx] <- pr_locus[idx] + pr_distinct * pr
    }

    # remove zeroes
    pr_locus <- pr_locus[pr_locus > 0]

    pr_by_locus[[locus]] <- pr_locus
  }

  # sum across loci to obtain result
  results = list(pf = compute_pr_of_sum(pr_by_locus))
  results$by_locus = pr_by_locus
  results$noa = as.numeric(names(results$pf))
  results$min = min(results$noa)
  results$max = max(results$noa)

  class(results) = "pf"

  return(results)
}
