#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector pr_sum(IntegerVector x1, NumericVector fx1,
                     IntegerVector x2, NumericVector fx2) {


  if (x1.size() < 1 || x2.size() < 1){
    Rcpp::stop("x1 and x2 need to have a length of at least 1");
  }
  if (x1.size() != fx1.size()){
    Rcpp::stop("length of x1 and fx1 are not equal");
  }
  if (x2.size() != fx2.size()){
    Rcpp::stop("length of x2 and fx2 are not equal");
  }

  long sum_min = x1[0] + x2[0];
  long sum_max = x1[x1.size() - 1] + x2[x2.size() - 1];

  NumericVector fx_sum(sum_max - sum_min + 1);
  CharacterVector x_sum(fx_sum.size());

  for (int i = 0; i < x_sum.size(); i++){
    x_sum[i] = std::to_string(sum_min + i);
  }
  fx_sum.attr("names") = x_sum;

  // determine non-zero indices for fx2
  std::vector<int> idx2;
  for (int i = 0; i < x2.size(); i++){
    if (fx2[i] > 0) idx2.push_back(i);
  }

  // convolution
  for (int i1 = 0; i1 < x1.size(); i1++){

    if (fx1[i1] == 0) continue;

    for (int i = 0; i < idx2.size(); i++){
      int i2 = idx2[i];

      long x = x1[i1] + x2[i2];
      double fx = fx1[i1] * fx2[i2];

      fx_sum[x - sum_min] += fx;
    }
  }

  return fx_sum;
}
