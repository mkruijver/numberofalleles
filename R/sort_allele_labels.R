sort_allele_labels <- function(unsorted_allele_labels){
  unsorted_allele_labels_numeric <- as.numeric(unsorted_allele_labels)
  sortable <- as.character(unsorted_allele_labels_numeric) == unsorted_allele_labels

  c(as.character(sort(as.numeric(unsorted_allele_labels[sortable]))), unsorted_allele_labels[!sortable])
}
