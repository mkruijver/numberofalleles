#include <Rcpp.h>
using namespace Rcpp;

// implements equation (6) of Tvedebrink, T., On the exact distribution
// of the number of alleles in DNA mixtures, International Journal of
// Legal Medicine (2014)
double S_recursive_internal(NumericVector &f, IntegerVector &alpha,
                            int n,
                            NumericVector &power_sums,
                            std::map<std::vector<int>, double> &memo){

  // check if we have computed this value before
  std::vector<int> alpha_slice(alpha.begin(), alpha.begin() + n);
  std::sort(alpha_slice.begin(), alpha_slice.end());

  auto it = memo.find(alpha_slice);
  if (it != memo.end()) {
    return it->second;
  }

  if (n == 1){

    double pr = 0.0;
    pr = power_sums[alpha[0]];
    memo[alpha_slice] = pr;

    return pr;
  }
  else{

    double term_one = power_sums[alpha[n - 1]];

    double term_two = S_recursive_internal(f, alpha, n - 1, power_sums, memo);

    double sum_terms = 0;

    for (int i = 0; i < n - 1; i++){
      alpha[i] += alpha[n - 1];

      double sum_term = S_recursive_internal(f, alpha, n - 1, power_sums, memo);
      sum_terms += sum_term;

      alpha[i] -= alpha[n - 1]; // restore state
    }

    double result = term_one * term_two - sum_terms;

    memo[alpha_slice] = result;

    return(result);
  }

}

// [[Rcpp::export]]
bool is_equal(NumericVector x, NumericVector y) {

  if (x.size() != y.size()) return false;

  return &x[0] == &y[0];
}

double S_recursive_hw(NumericVector f, IntegerVector alpha) {

  static NumericVector power_sums(0);

  // pre compute sums of powers for some speed up
  int max_power = sum(alpha) + 1;

  // initialise map for memoisation
  static std::map<std::vector<int>, double> memo;
  static NumericVector last_f = NumericVector::create();

  bool different_f = !is_equal(last_f, f);
  if (different_f){
    // Rcout << "different f!\n";

    memo.clear();
    last_f = f;
  }else{
    // Rcout << "same f!\n";

  }

  if (different_f || (power_sums.length() < (max_power + 1))){
    power_sums = NumericVector(max_power + 1);
    for (int k = 0; k <= max_power; k++){
      double sum = 0;

      for (int a = 0; a < f.size(); a++){
        sum += std::pow(f[a], k);
      }

      power_sums[k] = sum;
    }

  }

  return S_recursive_internal(f, alpha, alpha.size(), power_sums, memo);
}

IntegerVector ek_(IntegerVector& x, int k){
  IntegerVector y(x.size());

  y[k - 1] = 1;
  return y;
}

double S_recursive_fst_internal(NumericVector f, double fst,
                                IntegerVector a, IntegerVector b){

  // adapted from https://github.com/mikldk/DNAtools/blob/4abbd7a31de6407d411363a15cd72d22ac62c959/src/class_probsObj.h#L397

  int len_b = b.size();

  if(len_b == 0){
    return S_recursive_hw(f, a);
  }
  else if(is_true(any(b == 0))){
    return S_recursive_fst_internal(f, fst, a, b[b != 0]);
  }
  else{
    if(b[len_b - 1] == 1 && len_b == 1){
      return S_recursive_fst_internal(f, fst, a + ek_(a, len_b), head(b, len_b - 1));
    }
    else if(b[len_b - 1] == 1){
      // print(b);
      return (1.0 - fst) / (1 + (sum(b) - 2) * fst) *
        S_recursive_fst_internal(f, fst, a + ek_(a, len_b), head(b, len_b - 1));
    }
    else{
      return
       ((b[len_b - 1] - 1) * fst * S_recursive_fst_internal(f, fst, a, b - ek_(b, len_b)) +
        (1.0 - fst) * S_recursive_fst_internal(f, fst, a + ek_(a, len_b), b - ek_(b, len_b))) /
          (1 + (sum(b) - 2) * fst);
    }
  }
}

// [[Rcpp::export]]
double S_recursive(NumericVector f, IntegerVector alpha, double fst) {

  if (fst!=0 && fst < 1e-16){
    Rcpp::stop("fst!=0 && fst < 1e-16");
  }
  if (fst > 1){
    Rcpp::stop("fst > 1");
  }

  if (abs(fst) == 0){
    return S_recursive_hw(f, alpha);
  }
  else{
    IntegerVector a(alpha.size());

    return(S_recursive_fst_internal(f, fst, a, alpha));
  }
}
