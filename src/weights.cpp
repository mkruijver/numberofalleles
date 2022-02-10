#include <Rcpp.h>
using namespace Rcpp;

const unsigned long long factorials[] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800,
                             479001600, 6227020800, 87178291200, 1307674368000, 20922789888000,
                             355687428096000, 6402373705728000, 121645100408832000, 2432902008176640000};


unsigned long long get_factorial(int x){
  if (x < 0){
    Rcpp::stop("can not obtain factorial for negative number");
  }
  if (x > 20){
    Rcpp::stop("can not obtain factorial for x > 20");
  }

  return factorials[x];
}

// [[Rcpp::export]]
unsigned long long weights_cpp(IntegerVector x) {

  if (x.size() == 0) return 0;
  if (x.size() == 1) return 1;

  for (int i = 0; i < x.size(); i++){
    if (x[i]==0){
      Rcpp::stop("cannot compute weights for alpha with zeroes");
    }
  }

  unsigned long long sum_x = 0;
  unsigned long long prod_factorial_x = 1;
  for (int i = 0; i < x.size(); i++){
    sum_x += x[i];
    prod_factorial_x *= get_factorial(x[i]);
  }

  // compute prod(factorial(rle(sort(x))$lengths))
  unsigned long long prod_factorial_rle  = 1;

  std::sort(x.begin(), x.end());
  int run_length = 1;

  for (int i = 1; i < x.size(); i++){
    if (x[i-1] == x[i]){
      run_length++;
    }
    else{
      // Rcout << "run of " << run_length << "\n";
      prod_factorial_rle *= get_factorial(run_length);

      run_length = 1;
    }
  }
  prod_factorial_rle *= get_factorial(run_length);

  return get_factorial(sum_x) / (prod_factorial_x * prod_factorial_rle);
}
// get_weights <- function(x){
//   x <- x[x>0]
//
//   factorial(sum(x)) /
//     (prod(factorial(x)) * prod(factorial(rle(sort(x))$lengths)))
// }
