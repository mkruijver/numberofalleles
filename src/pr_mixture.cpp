#include <Rcpp.h>
using namespace Rcpp;

const double factorials[] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600,
                             6227020800,87178291200,1307674368000,20922789888000,355687428096000,6402373705728000,
                             121645100408832000, 2432902008176640000};

double pr_mixture_given_n_independent_alleles_numerator(std::vector<int> &ri, int &c, std::vector<double> &f, double &theta, IntegerVector &R){
  double numerator = 1.;

  for(int i = 0; i < c; ++i){
    for(int j = 0; j <= ri[i]; ++j){
      numerator *= ((1. - theta) * f[R[i]-1] + j * theta);
    }
  }

  return(numerator);
}

void pr_mixture_given_n_independent_alleles_recurs(std::vector<int> &ri, int j, int rleft, int &c,
               double &total_pr, double &theta, IntegerVector &R, std::vector<double> &f){
  if (j < (c-1)){
    for(int rj = 0; rj <= rleft; ++rj){
      ri[j] = rj;
      pr_mixture_given_n_independent_alleles_recurs(ri, j+1, rleft-rj, c, total_pr, theta, R, f);
    }
  }else{
    ri[j] = rleft;
    // compute the pr.
    double num = pr_mixture_given_n_independent_alleles_numerator(ri,c,f,theta,R);
    double denom=1.;

    for(int i = 0; i < c; ++i) denom *= factorials[ri[i]+1];

    total_pr += num/denom;
  }
}

// [[Rcpp::export]]
double pr_mixture_given_n_independent_alleles(int number_of_independent_alleles, double theta, IntegerVector R, std::vector<double> f){
  // computes the pr of observing a single replicate R ignoring the possibility of dropout / drop-in
  // This implements equation (1) in https://doi.org/10.1111/j.1556-4029.2010.01550.x

  int c = R.size(); // # alleles in R
  int r = number_of_independent_alleles - c;   // # `unconstrained' alleles

  if(number_of_independent_alleles > 20) Rcpp::stop("number_of_independent_alleles > 20 is not supported");
  if(number_of_independent_alleles < 1) Rcpp::stop("number_of_independent_alleles < 1 is not supported");

  // start recursion
  std::vector<int> ri(c);

  double pr = 0.;
  if (R.size()>0) pr_mixture_given_n_independent_alleles_recurs(ri, 0, r, c, pr, theta, R, f);

  // multiply with factor taken out of the sum
  double fact = factorials[number_of_independent_alleles];
  for(int j = 0; j < number_of_independent_alleles; ++j) fact/= ((1. - theta) + j * theta);

  return(pr * fact);
}

// [[Rcpp::export]]
NumericVector pr_mixtures_given_n_independent_alleles(int number_of_independent_alleles, double theta, IntegerMatrix mixtures, std::vector<double> f){
  int number_of_mixtures = mixtures.ncol();

  NumericVector prs(number_of_mixtures);

  for (int i_mixture = 0; i_mixture < number_of_mixtures; i_mixture++){

    IntegerMatrix::Column mix = mixtures( _ , i_mixture);

    prs[i_mixture] = pr_mixture_given_n_independent_alleles(number_of_independent_alleles, theta, mix, f);
  }

  return prs;
}

// [[Rcpp::export]]
NumericVector pr_k_allele_mixtures_given_n_independent_alleles(int number_of_independent_alleles, IntegerMatrix x, std::vector<double> f){
  // computes the pr of observing k-allele mixture ignoring the possibility of dropout / drop-in
  // if the k-allele mixture is sampled according to freqs without replacement

  NumericVector prs(x.ncol());

  int k = x.nrow();

  for (int i_mixture = 0; i_mixture < x.ncol(); i_mixture++){
    double pr = 1.0;

    double f_cum_pr = 1.0;

    for (int i_allele = 0; i_allele < k; i_allele++){
      int a = x(i_allele, i_mixture) - 1;
      pr *= (f[a] / f_cum_pr);

      f_cum_pr -= f[a];
    }

    prs[i_mixture] = pr * factorials[number_of_independent_alleles];
  }

  return prs;
}

