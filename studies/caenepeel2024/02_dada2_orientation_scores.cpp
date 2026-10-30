#include <Rcpp.h>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <vector>

using namespace Rcpp;

static int tax_kmer(const char *seq, unsigned int k) {
  int kmer = 0;
  for (unsigned int j = 0; j < k; ++j) {
    unsigned int nti;
    if (seq[j] == 'A') nti = 0;
    else if (seq[j] == 'C') nti = 1;
    else if (seq[j] == 'G') nti = 2;
    else if (seq[j] == 'T') nti = 3;
    else return -1;
    kmer = 4 * kmer + nti;
  }
  return kmer;
}

static std::vector<int> tax_karray(const std::string &seq, unsigned int k) {
  std::vector<int> out;
  if (seq.size() < k) return out;
  out.reserve(seq.size() - k + 1);
  for (size_t i = 0; i <= seq.size() - k; ++i) {
    int kmer = tax_kmer(seq.c_str() + i, k);
    if (kmer >= 0) out.push_back(kmer);
  }
  std::sort(out.begin(), out.end());
  return out;
}

static float best_log_score(
  const std::vector<int> &kmers,
  const std::vector<float> &log_probability,
  size_t n_kmers,
  size_t n_genus
) {
  float best = -FLT_MAX;
  for (size_t g = 0; g < n_genus; ++g) {
    const float *probability = &log_probability[g * n_kmers];
    float score = 0.0f;
    for (size_t pos = 0; pos < kmers.size(); ++pos) {
      score += probability[kmers[pos]];
      if (score < best) break;
    }
    if (best > 0.0f || score > best) best = score;
  }
  return best;
}

// [[Rcpp::export]]
DataFrame dada2_orientation_scores(
  CharacterVector sequences,
  CharacterVector reverse_complements,
  CharacterVector references,
  IntegerVector reference_to_genus
) {
  const unsigned int k = 8;
  const size_t n_kmers = 1u << (2 * k);
  const size_t n_sequences = sequences.size();
  const size_t n_references = references.size();

  if (reverse_complements.size() != n_sequences) {
    stop("sequences and reverse_complements must have equal lengths.");
  }
  if (reference_to_genus.size() != n_references) {
    stop("references and reference_to_genus must have equal lengths.");
  }
  if (n_references == 0) stop("At least one reference sequence is required.");

  int max_genus = max(reference_to_genus);
  if (max_genus < 1) stop("reference_to_genus must use positive one-based indices.");
  const size_t n_genus = static_cast<size_t>(max_genus);

  std::vector<float> genus_count_plus_one(n_genus, 1.0f);
  std::vector<float> kmer_prior(n_kmers, 0.0f);
  std::vector<float> log_probability(n_genus * n_kmers, 0.0f);
  std::vector<unsigned char> present(n_kmers, 0);

  for (size_t i = 0; i < n_references; ++i) {
    int g = reference_to_genus[i] - 1;
    if (g < 0 || static_cast<size_t>(g) >= n_genus) {
      stop("Invalid reference_to_genus index.");
    }
    genus_count_plus_one[g] += 1.0f;
    std::fill(present.begin(), present.end(), 0);
    std::string ref = as<std::string>(references[i]);
    if (ref.size() >= k) {
      for (size_t pos = 0; pos <= ref.size() - k; ++pos) {
        int kmer = tax_kmer(ref.c_str() + pos, k);
        if (kmer >= 0) present[kmer] = 1;
      }
    }
    float *genus_probability = &log_probability[static_cast<size_t>(g) * n_kmers];
    for (size_t kmer = 0; kmer < n_kmers; ++kmer) {
      if (present[kmer]) {
        genus_probability[kmer] += 1.0f;
        kmer_prior[kmer] += 1.0f;
      }
    }
    if ((i + 1) % 512 == 0) checkUserInterrupt();
  }

  for (size_t kmer = 0; kmer < n_kmers; ++kmer) {
    kmer_prior[kmer] = (kmer_prior[kmer] + 0.5f) /
      (1.0f + static_cast<float>(n_references));
  }
  for (size_t g = 0; g < n_genus; ++g) {
    float *genus_probability = &log_probability[g * n_kmers];
    for (size_t kmer = 0; kmer < n_kmers; ++kmer) {
      genus_probability[kmer] = std::log(
        (genus_probability[kmer] + kmer_prior[kmer]) /
        genus_count_plus_one[g]
      );
    }
  }

  NumericVector forward_score(n_sequences);
  NumericVector reverse_score(n_sequences);
  CharacterVector decision(n_sequences);
  for (size_t i = 0; i < n_sequences; ++i) {
    std::vector<int> forward_kmers = tax_karray(as<std::string>(sequences[i]), k);
    std::vector<int> reverse_kmers = tax_karray(as<std::string>(reverse_complements[i]), k);
    float forward = best_log_score(forward_kmers, log_probability, n_kmers, n_genus);
    float reverse = best_log_score(reverse_kmers, log_probability, n_kmers, n_genus);
    forward_score[i] = forward;
    reverse_score[i] = reverse;
    decision[i] = reverse > forward ? "reverse" : "forward";
    if ((i + 1) % 128 == 0) checkUserInterrupt();
  }

  return DataFrame::create(
    _["decision"] = decision,
    _["forward_log_score"] = forward_score,
    _["reverse_log_score"] = reverse_score,
    _["score_difference_reverse_minus_forward"] = reverse_score - forward_score,
    _["stringsAsFactors"] = false
  );
}
