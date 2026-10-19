#pragma once
#include "gsea/bayes_export.h"
#include "gsea/summary_enrichment.hpp"
#include <gbayes/bayesr.hpp>
namespace gsea {
struct BayesianPathwayOptions {
  std::string trait, response_definition, covariance_provenance;
  bool scale_response{true}, scale_features{true};
  std::size_t maximum_features{4096}, maximum_covariates{128}, maximum_block_size{2048};
  std::size_t maximum_preparation_bytes{512ULL*1024*1024};
  double minimum_reciprocal_condition{1e-10};
};
struct BayesianPathwayPreparation {
  std::size_t gene_count{}, nuisance_rank{}, block_count{};
  std::vector<std::string> gene_ids, nuisance_names;
  std::vector<double> feature_scales;
  double response_scale{1}, residual_sum_squares{}, nuisance_condition{}, worst_block_condition{};
  std::size_t estimated_preparation_bytes{}; // Matrix payload estimate, not process peak.
  BayesianPathwayOptions options;
};
struct PathwayPosterior {
  std::string name;
  double mean{}, variance{}, inclusion_probability{}; // Original feature/response units.
};
struct BayesianPathwayCResult {
  gbayes::BayesCResult native; // Sampler coordinates; preserves priors, chains and learning.
  BayesianPathwayPreparation preparation;
  std::vector<PathwayPosterior> pathways;
};
struct BayesianPathwayRResult {
  gbayes::BayesRResult native;
  BayesianPathwayPreparation preparation;
  std::vector<PathwayPosterior> pathways;
};
// Unlike set_design(), this builder imposes no frequentist coefficient-count limit.
GSEA_BAYES_EXPORT CompetitiveDesign bayesian_set_design(const PreparedSets&,
    const std::vector<std::string>& covariate_names={},const std::vector<double>& covariates={},
    const BayesianPathwayOptions& = {});
class GSEA_BAYES_EXPORT PreparedBayesianPathways {
public:
  // Binary sets with independent observations and an intercept only.
  // Uses direct set cross-products; preserves n-1 pathway normalization.
  PreparedBayesianPathways(const PreparedSets&, const NamedValues&,
      const BayesianPathwayOptions&);
  PreparedBayesianPathways(const CompetitiveDesign&, const NamedValues&,
      const std::vector<GeneCovarianceBlock>&, const BayesianPathwayOptions&);
  // Uses unadjusted transformed_response and adjusts BOTH response and features.
  PreparedBayesianPathways(const PreparedSummaryEnrichment&, const BayesianPathwayOptions&);
  const gbayes::FixedGaussianStatistics& statistics() const noexcept { return statistics_; }
  const gbayes::GaussianNormalization& normalization() const noexcept { return normalization_; }
  const BayesianPathwayPreparation& preparation() const noexcept { return preparation_; }
  BayesianPathwayCResult fit_bayesc(const gbayes::BayesCPrior& = {},
      const gbayes::BayesCControls& = {}, gbayes::GaussianLearning = {}) const;
  BayesianPathwayRResult fit_bayesr(const gbayes::BayesRPrior& = {},
      const gbayes::BayesCControls& = {}, gbayes::GaussianLearning = {}) const;
private:
  gbayes::FixedGaussianStatistics statistics_;
  gbayes::GaussianNormalization normalization_;
  BayesianPathwayPreparation preparation_;
};
}
