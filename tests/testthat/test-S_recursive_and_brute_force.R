test_that("S(f, alpha, fst) recursive and brute force are consistent", {

  set.seed(1)
  f0 <- runif(4)
  f <- f0 / sum(f0)

  for (fst in c(0,0.01,0.1,1)){

    for (m in 1:5){

      alphas <- apply(partitions::parts(m),2, function(x) x[x>0], simplify = FALSE)

      for (alpha in alphas){
        result_recursive <- numberofalleles:::S_recursive(f, alpha, fst)
        result_brute_force <- numberofalleles:::S_brute_force(f, alpha, fst)

        expect_equal(result_recursive, result_brute_force)
      }
    }
  }


})
