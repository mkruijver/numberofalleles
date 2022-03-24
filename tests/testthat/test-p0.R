test_that("p0 is consistent with HW calculation", {

  set.seed(123)
  f0 <- runif(10)
  f <- f0 / sum(f0)

  # HW
  p0 <- numberofalleles:::compute_p0(f, 5, fst = 0)
  expect_equal(p0, (1-f)^5)
})

test_that("p0 works for Fst=1", {

  set.seed(123)
  f0 <- runif(10)
  f <- f0 / sum(f0)

  expect_equal(numberofalleles:::compute_p0(f, 4, fst = 1), 1 - f)
})
