
compute_pr_of_sum <- function(x){
  if (length(x) == 1){
    return(x[[1]])
  }else if (length(x) == 2){

    return(pr_sum(x1 = as.integer(names(x[[1]])),
                  fx1 = x[[1]],
                  x2 = as.integer(names(x[[2]])),
                  fx2 = x[[2]]))

  }else{

    pr <- x[[1]]
    for (i in 2:length(x)){
      pr <- compute_pr_of_sum(list(pr, x[[i]]))
    }

    return(pr)
  }
}
