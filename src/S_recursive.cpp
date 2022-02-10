#include <Rcpp.h>
using namespace Rcpp;

// implements equation (6) of Tvedebrink, T., On the exact distribution
// of the number of alleles in DNA mixtures, International Journal of
// Legal Medicine (2014)
double S_recursive_internal(NumericVector &f, IntegerVector &alpha, int n){

  if (n == 1){

    double pr = 0.0;

    for (int i = 0; i < f.size(); i++){
      pr += std::pow(f[i], alpha[0]);
    }

    return pr;
  }
  else{

    double term_one = 0;
    for (int i = 0; i < f.size(); i++){
      term_one += std::pow(f[i], alpha[n - 1]);
    }

    double term_two = S_recursive_internal(f, alpha, n - 1);

    double sum_terms = 0;

    for (int i = 0; i < n - 1; i++){
      alpha[i] += alpha[n - 1];

      double sum_term = S_recursive_internal(f, alpha, n - 1);
      sum_terms += sum_term;

      alpha[i] -= alpha[n - 1]; // restore state
    }

    double result = term_one * term_two - sum_terms;

    return(result);
  }

}

// [[Rcpp::export]]
double S_recursive(NumericVector f, IntegerVector alpha) {
  return S_recursive_internal(f, alpha, alpha.size());
}
