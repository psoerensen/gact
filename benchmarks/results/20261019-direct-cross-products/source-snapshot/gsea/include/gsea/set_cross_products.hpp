#pragma once
#include "gsea/gsea.hpp"
namespace gsea {
// Binary sets, independent observations and optional intercept removal.
// No gene covariance, nuisance covariates or feature matrix is implicit.
struct SetCrossProducts {
  std::size_t n{};
  std::vector<std::string> feature_names;
  std::vector<double> xx, xy, centers, scales; // XX column-major
  double yy{}, response_center{}, response_scale{1};
};
GSEA_EXPORT SetCrossProducts set_cross_products(const PreparedSets&,
    const NamedValues&, bool center=true, bool scale_features=true,
    bool scale_response=true, std::size_t maximum_matrix_bytes=256ULL*1024*1024);
}
