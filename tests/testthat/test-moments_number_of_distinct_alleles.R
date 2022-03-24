test_that("expected value is consistent with by hand calculations", {

  set.seed(123)
  f0 <- runif(10)
  f <- f0 / sum(f0)

  for (m in 1:5){
    for (fst in c(0, 0.01, 0.3, 0.5, 1)){

      ev <- expected_number_of_distinct_alleles(
        number_of_independent_alleles = m, f = f, fst = fst)

      # via probability distribution
      p <- pr_number_of_distinct_alleles(number_of_independent_alleles = m,
                                         f = f, fst = fst)

      ev_manual <- as.vector(sum(as.numeric(names(p)) * p))

      expect_equal(ev, ev_manual)
    }
  }

})

test_that("variance is consistent with by hand calculations", {

  set.seed(2)
  f0 <- runif(10)
  f <- f0 / sum(f0)

  for (m in 1:5){
    for (fst in c(0, 0.01, 0.3, .99, 1)){

      var <- variance_number_of_distinct_alleles(
        number_of_independent_alleles = m, f = f, fst = fst)

      # via probability distribution
      p <- pr_number_of_distinct_alleles(number_of_independent_alleles = m,
                                         f = f, fst = fst)

      n <- as.numeric(names(p))
      ev <- as.vector(sum(n * p))

      var_manual <- sum(p * (n-ev)^2)

      expect_equal(var, var_manual)
    }
  }

})
