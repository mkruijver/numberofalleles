
test_that("var.pf is consistent with by hand calculations", {

  set.seed(3)
  f0 <- runif(10)
  f <- f0 / sum(f0)

  for (k in 1:5){
    for (fst in c(0, 0.01, 0.3, .99, 1)){

      p <- numberofalleles::pr_total_number_of_distinct_alleles(
        contributors = paste0("U", 1:k),
        freqs = list(locus1=f), fst = fst)

      var <- var(p)

      ev_manual <- as.vector(sum(p$noa * p$pf))
      var_manual <- sum(p$pf * (p$noa - ev_manual)^2)

      expect_equal(var, var_manual)
    }
  }

})
