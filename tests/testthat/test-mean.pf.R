test_that("mean.pf is consistent with by hand calculations", {

  set.seed(123)
  f0 <- runif(10)
  f <- f0 / sum(f0)

  for (k in 1:5){
    for (fst in c(0, 0.01, 0.3, 0.5, 1)){

      p <- numberofalleles::pr_total_number_of_distinct_alleles(
        contributors = paste0("U", 1:k),
        freqs = list(locus1=f), fst = fst)

      ev <- mean(p)

      # via probability distribution
      ev_manual <- as.vector(sum(p$noa * p$pf))

      expect_equal(ev, ev_manual)
    }
  }

})
