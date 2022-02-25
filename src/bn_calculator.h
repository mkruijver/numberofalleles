#include <Rcpp.h>
using namespace Rcpp;

class bn_calculator{

private:
  NumericVector _gamma0;
  double _gamma0_dot;

  NumericVector _gamma;
  double _gamma_dot;

public:
  int number_of_alleles;

  bn_calculator(NumericVector f, double fst){
    if (fst < 1e-16){
      Rcpp::stop("fst < 1e-16");
    }
    if (fst > 1){
      Rcpp::stop("fst > 1");
    }

    number_of_alleles = f.length();

    double alpha = (1.0 - fst) / fst;

    _gamma0 = alpha * f;
    _gamma0_dot = sum(_gamma0);

    reset();
  }

  double pr_next(int a){
    double pr = _gamma[a] / _gamma_dot;

    // increase pseudo count
    _gamma[a]++;
    _gamma_dot++;

    return pr;
  }

  void reset(){
    _gamma = clone(_gamma0);
    _gamma_dot = _gamma0_dot;
  }
};
