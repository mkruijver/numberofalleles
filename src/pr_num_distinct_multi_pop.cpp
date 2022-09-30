#include <Rcpp.h>
using namespace Rcpp;

const double factorials[] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600,
                             6227020800,87178291200,1307674368000,20922789888000,355687428096000,6402373705728000,
                             121645100408832000, 2432902008176640000};

void recurs_alleles(int i_allele, const int number_of_alleles, const IntegerVector &i_pop_by_allele,
                    const NumericMatrix &f, NumericVector &pr_num_distinct, uint64_t mask, double pr) {

  if (i_allele == number_of_alleles){
    // tail
    int number_distinct = __builtin_popcountll(mask);
    pr_num_distinct[number_distinct - 1] += pr;
  }
  else{
    double pr_backup = pr;
    uint64_t mask_backup = mask;

    for (int a = 0; a < f.nrow(); a++){
      double pr_a = f(a, i_pop_by_allele[i_allele]);

      mask |= (((uint64_t) 1) << a);
      pr *= pr_a;

      recurs_alleles(i_allele + 1, number_of_alleles, i_pop_by_allele, f, pr_num_distinct, mask, pr);

      mask = mask_backup;
      pr = pr_backup;
    }
  }
}

// [[Rcpp::export]]
NumericVector pr_num_distinct_multi_pop_brute_force(NumericMatrix f, IntegerVector num_indep_by_pop) {
  // this is a true brute force method that considers each possible allelci type for each indep allele

  int number_of_populations = f.ncol();

  if (num_indep_by_pop.size() != number_of_populations){
    Rcpp::stop("Length of num_indep_by_pop needs to equal number of columns of f");
  }

  int total_number_of_indep_alleles = sum(num_indep_by_pop);

  Function r_rep("rep");
  IntegerVector i_pop_by_allele = r_rep(Rcpp::seq(0, number_of_populations - 1), num_indep_by_pop);

  NumericVector pr_num_distinct(total_number_of_indep_alleles);
  recurs_alleles(0, total_number_of_indep_alleles, i_pop_by_allele, f, pr_num_distinct, 0, 1.0);

  pr_num_distinct.names() = Rcpp::seq(1, total_number_of_indep_alleles);

  return pr_num_distinct;
}

struct recurs_data_unordered_subsets{
  int number_of_alleles;
  int number_of_indep_alleles;

  double multinom_denom;

  uint64_t mask;
  double pr;

  double n_factorial;

  NumericVector f;
  std::vector<uint64_t> masks;
  std::vector<double> prs;
};

void recurs_and_create_unordered_subsets(int i_allele, int last_allele, int last_run_length, recurs_data_unordered_subsets &r){
  if (i_allele == r.number_of_indep_alleles){
    // tail
    r.masks.push_back(r.mask); // note that many of these masks are duplicates
    // duplicates could be removed by using a map but that is very expensive
    // it could be possible to generate subsets in lexicographical order instead
    r.prs.push_back(r.pr * (r.n_factorial / r.multinom_denom));

    // Rcpp::Rcout << "num distinct: " << __builtin_popcountll(r.mask) << " orders: " << r.multinom_denom << "\n";
  }
  else{
    uint64_t mask_backup = r.mask;
    double pr_backup = r.pr;
    double multinom_denom_backup = r.multinom_denom;

    for (int a = last_allele; a < r.number_of_alleles; a++){

      if (r.f[a] == 0) continue;

      r.pr *= r.f[a];
      r.mask |= (((uint64_t) 1) << a);

      if (a == last_allele){

        r.multinom_denom *= (last_run_length + 1);

        recurs_and_create_unordered_subsets(i_allele + 1, a, last_run_length + 1, r);
      }
      else{
        recurs_and_create_unordered_subsets(i_allele + 1, a, 1, r);
      }

      r.pr = pr_backup;
      r.mask = mask_backup;
      r.multinom_denom = multinom_denom_backup;
    }
  }
}

void recurs_cartesian_product(int i_population, uint64_t mask, double pr,
                              std::vector<recurs_data_unordered_subsets> &subsets,
                              NumericVector &pr_num_distinct){

  int number_of_populations = subsets.size();

  if (i_population == number_of_populations){
    // tail of recursion

    int number_distinct = __builtin_popcountll(mask);
    pr_num_distinct[number_distinct - 1] += pr;

  }
  else{
    uint64_t mask_backup = mask;
    double pr_backup = pr;

    std::vector<uint64_t> &masks = subsets[i_population].masks;
    std::vector<double> &prs = subsets[i_population].prs;

    for (int i_mask = 0; i_mask < masks.size(); i_mask++){
      mask |= masks[i_mask];
      pr *= prs[i_mask];

      recurs_cartesian_product(i_population + 1, mask, pr, subsets, pr_num_distinct);

      mask = mask_backup;
      pr = pr_backup;
    }

  }
}

// [[Rcpp::export]]
NumericVector pr_num_distinct_multi_pop(NumericMatrix f, IntegerVector num_indep_by_pop) {
  // this method is a bit more efficient because it precomputes masks
  // for all unordered allele combinations in each population
  // then takes the Cartesian product and determines the # distinct alleles for each

  int number_of_populations = f.ncol();

  if (num_indep_by_pop.size() != number_of_populations){
    Rcpp::stop("Length of num_indep_by_pop needs to equal number of columns of f");
  }

  int total_number_of_indep_alleles = sum(num_indep_by_pop);

  NumericVector pr_num_distinct(total_number_of_indep_alleles);
  pr_num_distinct.names() = Rcpp::seq(1, total_number_of_indep_alleles);

  int number_of_alleles = f.nrow();

  std::vector<recurs_data_unordered_subsets> subsets_by_pop(num_indep_by_pop.size());
  // List debug_subsets_by_pop(num_indep_by_pop.size());

  for (int i_pop = 0; i_pop < num_indep_by_pop.size(); i_pop++){
    int number_of_indep_alleles = num_indep_by_pop[i_pop];
    NumericVector f_pop = f(_ ,i_pop);

    recurs_data_unordered_subsets r{number_of_alleles, number_of_indep_alleles};
    r.f = f_pop;
    r.pr = 1.0;
    r.n_factorial = factorials[number_of_indep_alleles];
    r.multinom_denom = 1.0;

    recurs_and_create_unordered_subsets(0, 0, 0, r);

    subsets_by_pop[i_pop] = r;

    // debug_subsets_by_pop[i_pop] = List::create( _["pr"] = r.prs,  _["mask"] = r.masks);
  }

  // now recurs to sum over the Cartesian product of all subsets by pop
  recurs_cartesian_product(0, 0, 1.0, subsets_by_pop, pr_num_distinct);

  return pr_num_distinct;
}
