#' @title Compute the probability distribution of the number of independent alleles in a mixture with dropout and related contributors
#'
#' @param dropout_prs Numeric vector. Dropout probabilities per contributor.
#' @param ped_contributors Character vector with unique names of contributors. Valid names are the names of pedigree members.
#' @param pedigree \link[pedtools]{ped} object
#' @description When mixture contributors are related according to a pedigree, they may share some alleles identical by descent so that their total number of \emph{independent} alleles is smaller than two times the number of contributors. The number of \emph{independent} alleles can be further reduced if dropout plays a role. This function computes the probability distribution of the number of independent alleles that related mixture contributors have in total for a locus given their dropout parameters. Note that the number of \emph{distinct} alleles that is observed at the locus is typically smaller than the number of independent alleles due to allele sharing.
#' @seealso \link{pr_independent_alleles}
#' @examples
#' # without dropout, a father-mother-child mixture has 4 indep. alleles
#' p <- pr_independent_alleles_ped(pedtools::nuclearPed(),
#'                                 ped_contributors = as.character(1:3),
#'                                 dropout_prs = rep(0, 3))
#' stopifnot(identical(p,
#'                     stats::setNames(1., "4")))
#' @export
pr_independent_alleles_ped <- function(pedigree, ped_contributors,
                                       dropout_prs){

  if (!is.numeric(dropout_prs)){
    stop("dropout_prs needs to be a numeric vector")
  }

  if (any(dropout_prs < 0)){
    stop("dropout_prs needs to be non-negative")
  }

  if (any(dropout_prs > 1)){
    stop("dropout_prs should not exceed 1")
  }

  if (!inherits(pedigree, "ped")){
    stop("pedigree should be of class ped")
  }

  if (!all(ped_contributors %in% pedigree$ID)){
    stop("not all ped_contributors are found in pedigree")
  }

  if (length(ped_contributors) != length(dropout_prs)){
    stop("ped_contributors and dropout_prs do not have the same length")
  }

  if (!is.character(ped_contributors)){
    stop("ped_contributors is not a character vector")
  }

  multi_person_IBD <- ribd::multiPersonIBD(pedigree, ids = ped_contributors)

  # enumerate, for each person, the possible indep. alleles that are
  # dropped/not dropped out and their probs
  # there are 4 possibilities for each state: both drop, 1 drop (2x) and no drop
  drop_IBD_by_person <- setNames(
    sapply(seq_along(ped_contributors), function(i_contributor){
      prepare_drop_IBD(multi_person_IBD,
                       person_name = ped_contributors[i_contributor],
                       dropout_pr = dropout_prs[i_contributor])
    },
    simplify = FALSE), ped_contributors)

  # now iterate over the Cartesian product of all dropout probabilities
  # and sum the pr. of 0, 1, ... indep. alleles
  # for each multi person IBD state, there are 4 ^ (# contribs) possibilities
  ped_dropout <- data.frame(number_of_independent_alleles = 0:
                              (2*length(ped_contributors)))
  ped_dropout$pr <- 0. # we will add to this number in the loop

  ij <- do.call(expand.grid,
                replicate(n = length(ped_contributors), expr = 1:4, simplify = FALSE))
  ij_columns <- columns(ij)

  # loop over all multi person IBD states
  for (i_ibd in seq_len(nrow(multi_person_IBD))){

    # pr of dropout states for all contribs
    pr <- Reduce("*", lapply(seq_along(ped_contributors), function(i_contributor){
      drop_IBD_by_person[[i_contributor]][[i_ibd]][ij_columns[[i_contributor]],2]
    }))

    # corresponding number of unique alleles across contribs
    number_of_unique_alleles <- do.call(function(...)
      mapply(get_number_of_unique_alleles, ...),
      sapply(seq_along(ped_contributors), function(i_contributor){
        drop_IBD_by_person[[i_contributor]][[i_ibd]][ij_columns[[i_contributor]],3]}, simplify = FALSE))

    # collect sum of pr's by number of unique alleles
    ped_dropout$pr <- ped_dropout$pr +
      multi_person_IBD$Prob[i_ibd] *
      sapply(ped_dropout$number_of_independent_alleles, function(n){
        sum(pr[number_of_unique_alleles==n])
      })
  }

  ped_dropout <- ped_dropout[ped_dropout$pr > 0 ,]

  x <- ped_dropout$number_of_independent_alleles
  fx <- ped_dropout$pr

  stats::setNames(fx, x)
}

# this function finds, for the selected person, for each multi person IBD state,
# the possible sets of alleles that are left after some have dropped out
# with the corresponding probabilities

# example:
# ped <- pedtools::nuclearPed(father = "F", mother = "M",
#                     children = c("S1", "S2"))
# multi_person_ibd <- ribd::multiPersonIBD(ped, ids = c("S1", "S2"))
# prepare_drop_IBD(multi_person_ibd, "S2", 0.5)
prepare_drop_IBD <- function(multi_person_IBD, person_name, dropout_pr){

  drop_IBD_by_multi_person_IBD <- list()

  contrib_genos <- matrix(unlist(strsplit(multi_person_IBD[[person_name
  ]],
  split = " ", fixed = TRUE)),
  ncol = 2, byrow = TRUE)

  for(i_row in seq_len(nrow(multi_person_IBD))){
    pr <- multi_person_IBD[i_row, 1]

    a <- as.integer(contrib_genos[i_row, 1])
    b <- as.integer(contrib_genos[i_row, 2])

    dropout_prs <- c(dropout_pr^2,
                     dropout_pr * (1 - dropout_pr),
                     (1 - dropout_pr) * dropout_pr,
                     (1 - dropout_pr)^2)

    retained <- list(integer(), b, a, c(a,b))

    drop_IBD <- data.frame(geno = multi_person_IBD[[person_name
    ]][i_row],
    drop_pr = dropout_prs,
    retained = I(retained))

    names(drop_IBD) <- paste0(person_name, "_", names(drop_IBD))

    drop_IBD_by_multi_person_IBD[[1 + length(drop_IBD_by_multi_person_IBD)]] <-
      drop_IBD
  }

  drop_IBD_by_multi_person_IBD
}
