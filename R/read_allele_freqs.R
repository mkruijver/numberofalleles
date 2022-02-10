#' @title Read allele frequencies in FSIgen format (.csv)
#'
#' @param filename Path to csv file.
#' @param remove_zeroes Logical. Should frequencies of 0 be removed from the return value? Default is TRUE.
#' @param normalise Logical. Should frequencies be normalised to sum to 1? Default is TRUE.
#' @details Reads allele frequencies from a .csv file. The file should be in FSIgen format, i.e.
#'          comma separated with the first column specifying the allele labels and
#'          one column per locus. The last row should be the number of observations.
#'          No error checking is done since the file format is only loosely defined,
#'          e.g. we do not restrict the first column name or the last row name.
#' @return list
#' @examples
#' # below we read an allele freqs file that comes with the package
#' filename <- system.file("extdata","FBI_extended_Cauc.csv",package = "numberofalleles")
#' freqs <- read_allele_freqs(filename)
#' freqs # the output is just a list with an N attribute
#' @export
read_allele_freqs <- function(filename, remove_zeroes = TRUE, normalise = TRUE){
  raw <- readLines(filename)

  df_without_N <- utils::read.csv(file = filename,header = TRUE,
                         nrows = length(raw)-2, check.names=FALSE)
  df_with_N <- utils::read.csv(file = filename,
                      header = TRUE,nrows = length(raw)-1, check.names=FALSE)

  freqs <- list()
  alleles <- df_without_N[[1]]

  locus_idx <- seq(from=2,to=ncol(df_without_N),by=1)
  N <- stats::setNames(numeric(length(locus_idx)),nm = names(df_without_N[-1]))

  for(i_locus in locus_idx){
    f0 <- df_without_N[[i_locus]]
    f0[is.na(f0)] <- 0.
    f <- f0

    freqs[[names(df_without_N)[i_locus]]]  <- stats::setNames(f,nm =  alleles)
    N[names(df_without_N)[i_locus]] <- df_with_N[[i_locus]][length(df_with_N[[i_locus]])]
  }
  attr(freqs,"N") <- N

  if (remove_zeroes){
    freqs <- lapply(freqs, function(x) x[x>0])
  }

  if (normalise){
    freqs <- lapply(freqs, function(x) x/sum(x))
  }

  freqs
}
