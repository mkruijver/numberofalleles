#include <Rcpp.h>
using namespace Rcpp;

Rcpp::IntegerVector popcount(std::vector<uint64_t> x){

  int n = x.size();
  Rcpp::IntegerVector c(n);

  for (int i = 0; i < n; i++){
    c[i] = __builtin_popcountll(x[i]);
  }

  return c;
}

Rcpp::IntegerMatrix enumerate_genotypes(int number_of_alleles){

  int number_of_genotypes = number_of_alleles * (number_of_alleles + 1) / 2;

  Rcpp::IntegerMatrix g(number_of_genotypes, 2);

  int i_geno = 0;
  for (int a = 0; a < number_of_alleles; a++){
    for (int b = a; b < number_of_alleles; b++){
      g(i_geno, 0) = a;
      g(i_geno, 1) = b;
      i_geno++;
    }
  }

  return g;
}

std::vector<uint64_t> get_genotype_masks(Rcpp::IntegerMatrix g){

  if (g.ncol()!=2) Rcpp::stop("g needs to have 2 columns");

  int number_of_genotypes = g.nrow();

  std::vector<uint64_t> masks(number_of_genotypes);

  for (int i_geno = 0; i_geno < number_of_genotypes; i_geno++){
    uint64_t mask = 0;

    mask = ((uint64_t) 1) << g(i_geno, 0);
    mask |= ((uint64_t) 1) << g(i_geno, 1);

    masks[i_geno] = mask;
  }

  return masks;
}

std::vector<uint64_t> get_allele_masks(int number_of_alleles){
  std::vector<uint64_t> masks(number_of_alleles);

  for (int i_allele = 0; i_allele < number_of_alleles; i_allele++){
    masks[i_allele] = ((uint64_t) 1) << i_allele;
  }

  return masks;
}

Rcpp::NumericVector get_genotype_probabilities(Rcpp::IntegerMatrix g, Rcpp::NumericVector f){

  if (g.ncol()!=2) Rcpp::stop("g needs to have 2 columns");

  int number_of_genotypes = g.nrow();

  Rcpp::NumericVector genotype_probabilities(number_of_genotypes);

  for (int i_geno = 0; i_geno < number_of_genotypes; i_geno++){
    int a = g(i_geno, 0);
    int b = g(i_geno, 1);

    if (a==b){
      genotype_probabilities[i_geno] = f[a] * f[a];
    }
    else{
      genotype_probabilities[i_geno] = 2 * f[a] * f[b];
    }
  }

  return genotype_probabilities;
}

Rcpp::IntegerVector get_set_bits(uint64_t mask){
  int number_of_set_bits = __builtin_popcountll(mask);

  Rcpp::IntegerVector bits(number_of_set_bits);

  int i_bit = 0;

  while (mask != 0) {
    uint64_t t = mask & -mask;

    int r = __builtin_ctzll(mask);
    bits[i_bit++] = r;

    mask ^= t;
  }

  return bits;
}

Rcpp::IntegerVector get_set_bits_plus_one(uint64_t mask){
  int number_of_set_bits = __builtin_popcountll(mask);

  Rcpp::IntegerVector bits(number_of_set_bits);

  int i_bit = 0;

  while (mask != 0) {
    uint64_t t = mask & -mask;

    int r = __builtin_ctzll(mask);
    bits[i_bit++] = r + 1;

    mask ^= t;
  }

  return bits;
}

uint64_t get_mask(Rcpp::IntegerVector set_bits){
  uint64_t mask = 0;

  for (int i_bit = 0; i_bit < set_bits.size(); i_bit++){
    mask |= ((uint64_t) 1) << (set_bits[i_bit]);
  }

  return mask;
}

Rcpp::List masks_to_mixtures(std::vector<uint64_t> &masks, CharacterVector &allele_labels){
  int n = masks.size();

  Rcpp::List mixtures(n);

  for (int i = 0; i < n; i++){
    Rcpp::IntegerVector bits = get_set_bits_plus_one(masks[i]);

    bits.attr("class") = "factor";
    bits.attr("levels") = allele_labels;

    mixtures[i] = bits;
  }

  return mixtures;
}

IntegerVector r_sample_int(NumericVector pr, int number_of_samples, bool replace){
  Function sample_int("sample.int");

  return sample_int(_["n"] = pr.size(),
                    _["size"] = number_of_samples,
                    _["replace"] = replace,
                    _["prob"] = pr);
}

IntegerVector r_sample_int_uniform(int n, int number_of_samples, bool replace){
  Function sample_int("sample.int");

  return sample_int(_["n"] = n,
                    _["size"] = number_of_samples,
                    _["replace"] = replace);
}

std::vector<uint64_t> sample_mixture_masks_locus(NumericVector freqs_locus, int number_of_contributors, int number_of_samples) {

  IntegerMatrix g = enumerate_genotypes(freqs_locus.size());

  std::vector<uint64_t> masks = get_genotype_masks(g);

  NumericVector g_pr = get_genotype_probabilities(g, freqs_locus);

  std::vector<uint64_t> mixture_masks(number_of_samples);

  for (int i_contributor = 0; i_contributor < number_of_contributors; i_contributor++){
    IntegerVector s = r_sample_int(g_pr, number_of_samples, true);

    for (int i_sample = 0; i_sample < number_of_samples; i_sample++){
      mixture_masks[i_sample] |= masks[s[i_sample] - 1];
    }
  }

  return mixture_masks;
}

// [[Rcpp::export]]
std::vector<uint64_t> sample_mixture_masks_number_of_independent_alleles_locus(
    NumericVector freqs_locus, int number_of_indep_alleles, int number_of_samples) {

  std::vector<uint64_t> masks = get_allele_masks(freqs_locus.size());

  std::vector<uint64_t> mixture_masks(number_of_samples);

  for (int i_allele = 0; i_allele < number_of_indep_alleles; i_allele++){
    IntegerVector s = r_sample_int(freqs_locus, number_of_samples, true);

    for (int i_sample = 0; i_sample < number_of_samples; i_sample++){
      mixture_masks[i_sample] |= masks[s[i_sample] - 1];
    }
  }

  return mixture_masks;
}

// [[Rcpp::export]]
List sample_mixture_alleles_locus(NumericVector freqs_locus, int number_of_contributors, int number_of_samples) {

  std::vector<uint64_t> mixture_masks = sample_mixture_masks_locus(freqs_locus, number_of_contributors, number_of_samples);

  CharacterVector allele_labels = freqs_locus.attr("names");

  return masks_to_mixtures(mixture_masks, allele_labels);
}

// [[Rcpp::export]]
IntegerVector sample_allele_count_locus(NumericVector freqs_locus, int number_of_contributors, int number_of_samples) {

  std::vector<uint64_t> mixture_masks = sample_mixture_masks_locus(freqs_locus, number_of_contributors, number_of_samples);

  return popcount(mixture_masks);
}

// [[Rcpp::export]]
IntegerVector sample_allele_count_num_indep_alleles_locus(NumericVector freqs_locus, int number_of_indep_alleles, int number_of_samples) {

  std::vector<uint64_t> mixture_masks = sample_mixture_masks_number_of_independent_alleles_locus(freqs_locus, number_of_indep_alleles, number_of_samples);

  return popcount(mixture_masks);
}

// [[Rcpp::export]]
List sample_k_allele_mixtures_brute_force(NumericVector freqs_locus, int number_of_alleles, int number_of_samples){
  CharacterVector allele_labels = freqs_locus.attr("names");

  Function sample_int("sample.int");

  List mixtures(number_of_samples);

  for (int i_sample = 0; i_sample < number_of_samples; i_sample++){
    IntegerVector s = r_sample_int(freqs_locus, number_of_alleles, false);

    s.attr("class") = "factor";
    s.attr("levels") = allele_labels;

    mixtures[i_sample] = s;
  }

  return mixtures;
}

// [[Rcpp::export]]
IntegerMatrix sample_k_allele_mixtures(NumericVector freqs_locus, int number_of_alleles, int number_of_samples){

  // this is the first batch of alleles used in a rejection sampling scheme
  IntegerVector sampled_alleles = r_sample_int(freqs_locus, number_of_samples, true);

  IntegerMatrix samples(number_of_alleles, number_of_samples);

  int i_sampled_allele = 0;
  for (int i_sample = 0; i_sample < number_of_samples; i_sample++){

    uint64_t sample_mask = 0;

    int number_sampled = 0;

    while (number_sampled != number_of_alleles){
      if (i_sampled_allele == sampled_alleles.size()){
        // need a new batch
        sampled_alleles = r_sample_int(freqs_locus, number_of_samples, true);

        i_sampled_allele = 0;
      }

      // pop a sampled allele out of the batch
      int allele = sampled_alleles[i_sampled_allele] - 1;
      i_sampled_allele++;

      uint64_t allele_mask = ((uint64_t) 1) << allele;

      if ((sample_mask & allele_mask) == 0){
        samples(number_sampled, i_sample) = allele + 1;
        number_sampled++;
      }


            sample_mask |= allele_mask;
    }

  }

  return samples;
}

// [[Rcpp::export]]
IntegerMatrix sample_k_allele_mixtures_uniform(NumericVector freqs_locus, int number_of_alleles, int number_of_samples){

  // this is the first batch of alleles used in a rejection sampling scheme
  IntegerVector sampled_alleles = r_sample_int_uniform(freqs_locus.size(), number_of_samples, true);

  IntegerMatrix samples(number_of_alleles, number_of_samples);

  int i_sampled_allele = 0;
  for (int i_sample = 0; i_sample < number_of_samples; i_sample++){

    uint64_t sample_mask = 0;

    int number_sampled = 0;

    while (number_sampled != number_of_alleles){
      if (i_sampled_allele == sampled_alleles.size()){
        // need a new batch
        sampled_alleles = sampled_alleles = r_sample_int_uniform(freqs_locus.size(), number_of_samples, true);

        i_sampled_allele = 0;
      }

      // pop a sampled allele out of the batch
      int allele = sampled_alleles[i_sampled_allele] - 1;
      i_sampled_allele++;

      uint64_t allele_mask = ((uint64_t) 1) << allele;

      if ((sample_mask & allele_mask) == 0){
        samples(number_sampled, i_sample) = allele + 1;
        number_sampled++;
      }

      sample_mask |= allele_mask;
    }

  }

  return samples;
}
