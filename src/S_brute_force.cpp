#include "bn_calculator.h"

#include <Rcpp.h>
using namespace Rcpp;

void recurs_S_brute_force(int i_recurs, IntegerVector i,
                          NumericVector f, IntegerVector alpha,
                          std::vector<bool> &a_available,
                          double &total_pr){
  Rcpp::checkUserInterrupt();

  if (i_recurs == alpha.length()){
    // tail of recursion
    double pr = 1.0;

    for (int k = 0; k < i.length(); k++){
      pr *= std::pow(f[i[k]], alpha[k]);
    }

    total_pr += pr;
  }
  else{

    for (i[i_recurs] = 0; i[i_recurs] < f.length(); i[i_recurs]++){
      if (a_available[i[i_recurs]]){
        a_available[i[i_recurs]] = false;

        recurs_S_brute_force(i_recurs + 1, i, f, alpha, a_available, total_pr);

        a_available[i[i_recurs]] = true;
      }

    }
  }
}

double S_brute_force_hw(NumericVector f, IntegerVector alpha) {

  IntegerVector i(alpha.length());

  double total_pr = 0.0;

  std::vector<bool> a_available(f.size(), true);

  recurs_S_brute_force(0, i, f, alpha, a_available, total_pr);

  return total_pr;
}

void recurs_S_brute_force_fst(int i_recurs, IntegerVector i,
                              bn_calculator &bn, IntegerVector alpha,
                              std::vector<bool> &a_available,
                              double &total_pr){
  Rcpp::checkUserInterrupt();

  if (i_recurs == alpha.length()){
    // tail of recursion
    double pr = 1.0;

    // reset the BN counts
    bn.reset();

    for (int k = 0; k < i.length(); k++){
      for (int a = 0; a < alpha[k]; a++){

        pr *= bn.pr_next(i[k]);
      }
    }

    total_pr += pr;
  }
  else{

    for (i[i_recurs] = 0; i[i_recurs] < bn.number_of_alleles; i[i_recurs]++){
      if (a_available[i[i_recurs]]){
        a_available[i[i_recurs]] = false;

        recurs_S_brute_force_fst(i_recurs + 1, i, bn, alpha, a_available, total_pr);

        a_available[i[i_recurs]] = true;
      }

    }
  }
}

// [[Rcpp::export]]
double S_brute_force(NumericVector f, IntegerVector alpha, double fst = 0) {

  if (fst == 0){
    return S_brute_force_hw(f, alpha);
  }

  IntegerVector i(alpha.length());

  double total_pr = 0.0;

  std::vector<bool> a_available(f.size(), true);

  bn_calculator bn(f, fst);

  recurs_S_brute_force_fst(0, i, bn, alpha, a_available, total_pr);

  return total_pr;
}
