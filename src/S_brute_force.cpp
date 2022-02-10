#include <Rcpp.h>
using namespace Rcpp;

void recurs_S_brute_force(int i_recurs, IntegerVector i,
                          NumericVector f, IntegerVector alpha,
                          std::vector<bool> &a_available,
                          double &total_pr){

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

// [[Rcpp::export]]
double Sbruteforce(NumericVector f, IntegerVector alpha) {

  IntegerVector i(alpha.length());

  double total_pr = 0.0;

  std::vector<bool> a_available(f.size(), true);

  recurs_S_brute_force(0, i, f, alpha, a_available, total_pr);

  return total_pr;
}
