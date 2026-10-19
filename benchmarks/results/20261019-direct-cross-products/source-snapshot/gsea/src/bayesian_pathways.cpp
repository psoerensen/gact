#include "gsea/multitrait_pathways.hpp"
#include "gsea/set_cross_products.hpp"
#include <gsolve/dense_spd_factorization.hpp>
#include <gsolve/symmetric_eigensystem.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_set>
namespace gsea {
namespace {
void require(bool ok,const char* message) { if(!ok)throw std::invalid_argument(message); }
void names(const std::vector<std::string>& v) {
  std::unordered_set<std::string> seen;
  for(const auto& s:v)require(!s.empty()&&seen.insert(s).second,"Names must be nonempty and unique");
}
std::size_t product(std::size_t a,std::size_t b) {
  require(b==0||a<=std::numeric_limits<std::size_t>::max()/b,"Matrix extent overflow");return a*b;
}
void budget_add(std::size_t& bytes,std::size_t count,const BayesianPathwayOptions& o) {
  const auto extra=product(count,sizeof(double));
  require(bytes<=o.maximum_preparation_bytes&&extra<=o.maximum_preparation_bytes-bytes,"Bayesian preparation exceeds matrix payload budget");bytes+=extra;
}
std::vector<double> solve(const gsolve::DenseSpdFactorization& f,const std::vector<double>& b,std::size_t columns) {
  std::vector<double> out(b.size());const auto n=f.dimension();
  const auto status=f.solve({b.data(),n,columns,n},{out.data(),n,columns,n});
  if(!status)throw std::runtime_error(status.message());
  return out;
}
double condition(const std::vector<double>& a,std::size_t n,double minimum) {
  gsolve::SymmetricEigensystemOptions o;o.maximum_dimension=n;
  const auto e=gsolve::compute_symmetric_eigensystem({a.data(),n,n,n},o);
  require(e.status.is_ok(),"Invalid symmetric covariance/information matrix");
  const auto hi=e.eigenvalues.front(),lo=e.eigenvalues.back();
  require(hi>0&&lo>0&&lo/hi>=minimum,"Singular or ill-conditioned covariance/nuisance design; no automatic repair");return hi/lo;
}
CompetitiveDesign summary_design(const PreparedSummaryEnrichment& p) {
  auto d=p.design;
  if(p.options.residualize_internal&&!p.nuisance_design.feature_names.empty()) {
    require(p.nuisance_design.gene_ids==d.gene_ids,"Summary nuisance gene order differs");
    d.covariate_names.insert(d.covariate_names.end(),p.nuisance_design.feature_names.begin(),p.nuisance_design.feature_names.end());
    d.covariates.insert(d.covariates.end(),p.nuisance_design.features.begin(),p.nuisance_design.features.end());
  }
  return d;
}
gbayes::GaussianLearning learning_for(gbayes::GaussianLearning l,const gbayes::GaussianNormalization& n) {
  const auto& supplied=l.normalization;
  const bool empty=supplied.phenotype_sum_squares==0&&supplied.likelihood_dimension==0&&supplied.space==gbayes::GaussianSpace::full;
  require(empty||(supplied.phenotype_sum_squares==n.phenotype_sum_squares&&supplied.likelihood_dimension==n.likelihood_dimension&&supplied.space==n.space),"Normalization conflicts with prepared pathway likelihood");
  l.normalization=n;return l;
}
template<class Fit> std::vector<PathwayPosterior> posterior(const Fit& f,const BayesianPathwayPreparation& p) {
  std::vector<PathwayPosterior> out;out.reserve(f.marker_ids.size());
  for(std::size_t j=0;j<f.marker_ids.size();++j) {
    const double scale=p.response_scale/p.feature_scales[j];
    const double mean=f.mean[j]*scale,var=f.variance[j]*scale*scale;
    require(std::isfinite(mean)&&std::isfinite(var),"Posterior outside original-unit numerical range");
    out.push_back({f.marker_ids[j],mean,var,f.pip[j]});
  }
  return out;
}
}
CompetitiveDesign bayesian_set_design(const PreparedSets& p,const std::vector<std::string>& cn,
    const std::vector<double>& c,const BayesianPathwayOptions& o) {
  const auto n=p.universe().ids.size(),f=p.names().size();
  require(f>0&&f<=o.maximum_features&&cn.size()<=o.maximum_covariates,"Too many pathways/covariates");
  require(c.size()==product(n,cn.size()),"Covariate dimensions differ");
  std::size_t bytes=0;budget_add(bytes,product(n,f),o);
  CompetitiveDesign d{p.universe().ids,cn,p.names(),c,std::vector<double>(n*f,0)};
  for(std::size_t j=0;j<f;++j)for(auto k=p.offsets()[j];k<p.offsets()[j+1];++k)d.features[p.members()[k]+j*n]=1;
  return d;
}
PreparedBayesianPathways::PreparedBayesianPathways(const PreparedSummaryEnrichment& p,const BayesianPathwayOptions& o)
  :PreparedBayesianPathways(summary_design(p),p.transformed_response,p.covariance,o) {}
namespace {
struct ProjectedPathways {
  std::vector<gbayes::FixedGaussianStatistics> statistics;
  std::vector<BayesianPathwayPreparation> preparations;
  std::vector<double> response_gram;
  double dimension{};
};
// Shared scalar/multitrait preparation: one set of block factors, one nuisance
// projection and one pathway cross-product object, regardless of trait count.
ProjectedPathways project_pathways(const CompetitiveDesign& d,const std::vector<std::string>& response_ids,
    const std::vector<double>& response_values,
    const std::vector<std::string>& traits,const std::vector<GeneCovarianceBlock>& blocks,
    const BayesianPathwayOptions& o) {
  BayesianPathwayPreparation preparation;
  const auto t=traits.size();
  require(t>0&&t<=64,"Expected 1..64 traits in common-design preparation");
  names(traits);
  const auto n=d.gene_ids.size(),f=d.feature_names.size(),q=d.covariate_names.size()+1;
  require(f>0&&f<=o.maximum_features&&d.covariate_names.size()<=o.maximum_covariates&&n>q,"Invalid pathway/nuisance dimensions");
  require(d.features.size()==product(n,f)&&d.covariates.size()==product(n,q-1),"Design dimensions differ");
  require(response_ids==d.gene_ids&&response_values.size()==product(n,t),"Response gene identities/order differ");
  require(!o.response_definition.empty()&&!o.covariance_provenance.empty(),"Response definition and covariance provenance required");
  require(std::isfinite(o.minimum_reciprocal_condition)&&o.minimum_reciprocal_condition>0&&o.minimum_reciprocal_condition<1,"Invalid condition threshold");
  names(d.gene_ids);std::vector<std::string> all={"(Intercept)"};all.insert(all.end(),d.covariate_names.begin(),d.covariate_names.end());
  preparation.nuisance_names=all;all.insert(all.end(),d.feature_names.begin(),d.feature_names.end());names(all);
  for(const auto* v:{&d.features,&d.covariates,&response_values})for(double x:*v)require(std::isfinite(x),"Nonfinite response/design");
  require(!blocks.empty(),"Explicit covariance blocks required");
  std::size_t bytes=0,maxblock=0;std::vector<bool> seen(n,false);
  for(const auto& b:blocks) {
    const auto k=b.indices.size();require(k>0&&k<=o.maximum_block_size&&b.covariance.size()==product(k,k),"Invalid covariance block dimensions");
    for(auto i:b.indices){require(i<n&&!seen[i],"Covariance must partition every gene exactly once");seen[i]=true;}
    maxblock=std::max(maxblock,k);budget_add(bytes,product(4,product(k,k)),o);
  }
  require(std::all_of(seen.begin(),seen.end(),[](bool b){return b;}),"Covariance omits genes");
  // Conservative matrix-payload estimate, including factor/eigensolver scratch.
  budget_add(bytes,product(3,product(n,q+f+t)),o);
  budget_add(bytes,product(4,product(q+f+t,q+f+t)),o);
  budget_add(bytes,product(4,product(maxblock,q+f+t)),o);
  preparation.estimated_preparation_bytes=bytes;preparation.options=o;
  preparation.gene_ids=d.gene_ids;preparation.gene_count=n;preparation.nuisance_rank=q;preparation.block_count=blocks.size();
  std::vector<double> z(n,1);z.insert(z.end(),d.covariates.begin(),d.covariates.end());
  // Unit scaling improves nuisance conditioning without changing its column space.
  for(std::size_t j=1;j<q;++j){double s=0;for(std::size_t i=0;i<n;++i)s=std::max(s,std::abs(z[i+j*n]));require(s>0,"Zero nuisance column");for(std::size_t i=0;i<n;++i)z[i+j*n]/=s;}
  std::vector<double> r=response_values;r.insert(r.end(),d.features.begin(),d.features.end());
  // Scale first to protect products; record original units for final back-transform.
  std::vector<double> initial(f+t,0);
  for(std::size_t j=0;j<f+t;++j){for(std::size_t i=0;i<n;++i)initial[j]=std::max(initial[j],std::abs(r[i+j*n]));require(initial[j]>0,"Zero response/pathway column");for(std::size_t i=0;i<n;++i)r[i+j*n]/=initial[j];}
  std::vector<gsolve::DenseSpdFactorization> factors;factors.reserve(blocks.size());
  std::vector<double> zz(q*q,0),zr(q*(f+t),0);
  for(const auto& b:blocks) {
    const auto k=b.indices.size();preparation.worst_block_condition=std::max(preparation.worst_block_condition,condition(b.covariance,k,o.minimum_reciprocal_condition));
    factors.emplace_back(gsolve::ConstMatrixView{b.covariance.data(),k,k,k});require(factors.back().status().is_ok(),"Covariance factorization failed");
    std::vector<double> zb(k*q),rb(k*(f+t));
    for(std::size_t j=0;j<q;++j)for(std::size_t i=0;i<k;++i)zb[i+j*k]=z[b.indices[i]+j*n];
    for(std::size_t j=0;j<f+t;++j)for(std::size_t i=0;i<k;++i)rb[i+j*k]=r[b.indices[i]+j*n];
    const auto wz=solve(factors.back(),zb,q),wr=solve(factors.back(),rb,f+t);
    for(std::size_t j=0;j<q;++j)for(std::size_t i=0;i<=j;++i){long double v=0;for(std::size_t h=0;h<k;++h)v+=static_cast<long double>(zb[h+i*k])*wz[h+j*k];zz[i+j*q]+=static_cast<double>(v);}
    for(std::size_t j=0;j<f+t;++j)for(std::size_t i=0;i<q;++i){long double v=0;for(std::size_t h=0;h<k;++h)v+=static_cast<long double>(zb[h+i*k])*wr[h+j*k];zr[i+j*q]+=static_cast<double>(v);}
  }
  for(std::size_t j=0;j<q;++j)for(std::size_t i=0;i<j;++i)zz[j+i*q]=zz[i+j*q];
  preparation.nuisance_condition=condition(zz,q,o.minimum_reciprocal_condition);
  const gsolve::DenseSpdFactorization nuisance({zz.data(),q,q,q});require(nuisance.status().is_ok(),"Nuisance factorization failed");
  const auto coefficients=solve(nuisance,zr,f+t);
  for(std::size_t j=0;j<f+t;++j)for(std::size_t i=0;i<n;++i){long double v=r[i+j*n];for(std::size_t h=0;h<q;++h)v-=static_cast<long double>(z[i+h*n])*coefficients[h+j*q];r[i+j*n]=static_cast<double>(v);}
  std::vector<double> gram((f+t)*(f+t),0);
  for(std::size_t b=0;b<blocks.size();++b) {
    const auto& ix=blocks[b].indices;const auto k=ix.size();std::vector<double> rb(k*(f+t));
    for(std::size_t j=0;j<f+t;++j)for(std::size_t i=0;i<k;++i)rb[i+j*k]=r[ix[i]+j*n];
    const auto wr=solve(factors[b],rb,f+t);
    for(std::size_t j=0;j<f+t;++j)for(std::size_t i=0;i<=j;++i){long double v=0;for(std::size_t h=0;h<k;++h)v+=static_cast<long double>(rb[h+i*k])*wr[h+j*k];gram[i+j*(f+t)]+=static_cast<double>(v);}
  }
  const double dimension=static_cast<double>(n-q);std::vector<double> scale(f+t,1);
  for(std::size_t j=0;j<f+t;++j) {
    const double v=gram[j+j*(f+t)];
    if(!std::isfinite(v)||v<=1e-24*static_cast<double>(n))throw std::invalid_argument("No usable residual variation in "+(j<t?traits[j]:d.feature_names[j-t]));
    scale[j]=(j<t?o.scale_response:o.scale_features)?std::sqrt(v/dimension):1/initial[j];
    require(std::isfinite(scale[j])&&scale[j]>0,"Scaling outside numerical range");
  }

  for(std::size_t j=t;j<f+t;++j)preparation.feature_scales.push_back(initial[j]*scale[j]);
  for(double x:preparation.feature_scales)require(std::isfinite(x)&&x>0,"Feature scale overflow");
  std::vector<double> c(f*f);
  for(std::size_t j=0;j<f;++j)for(std::size_t i=0;i<=j;++i)
    c[i+j*f]=c[j+i*f]=gram[i+t+(j+t)*(f+t)]/scale[i+t]/scale[j+t];
  const auto common=gbayes::dense_cross_products(f,std::move(c));
  ProjectedPathways out;out.dimension=dimension;out.response_gram.resize(t*t);
  for(std::size_t j=0;j<t;++j)for(std::size_t i=0;i<=j;++i) {
    const auto value=gram[i+j*(f+t)]/scale[i]/scale[j];
    require(std::isfinite(value),"Normalized response cross-product overflow");
    out.response_gram[i+j*t]=out.response_gram[j+i*t]=value;
  }
  for(std::size_t h=0;h<t;++h) {
    preparation.options.trait=traits[h];
    preparation.response_scale=initial[h]*scale[h];
    preparation.residual_sum_squares=out.response_gram[h+h*t];
    require(std::isfinite(preparation.response_scale)&&preparation.response_scale>0,"Response scale overflow");
    require(preparation.residual_sum_squares>0,"Normalized response outside numerical range");
    gbayes::FixedGaussianStatistics statistics;statistics.marker_ids=d.feature_names;statistics.trait=traits[h];
    statistics.effect_units="prepared gene-response / prepared pathway feature";
    statistics.provenance=o.response_definition+"; "+o.covariance_provenance+"; GLS nuisance projection of response AND pathways";
    statistics.cross_products=common;statistics.score.resize(f);
    for(std::size_t j=0;j<f;++j) {
      statistics.score[j]=gram[h+(j+t)*(f+t)]/scale[h]/scale[j+t];
      require(std::isfinite(statistics.score[j]),"Normalized score overflow");
    }
    out.preparations.push_back(preparation);out.statistics.push_back(std::move(statistics));
  }
  return out;
}
}
PreparedBayesianPathways::PreparedBayesianPathways(const CompetitiveDesign& d,const NamedValues& y,
    const std::vector<GeneCovarianceBlock>& blocks,const BayesianPathwayOptions& o) {
  auto projected=project_pathways(d,y.ids,y.values,{o.trait},blocks,o);
  statistics_=std::move(projected.statistics.front());preparation_=std::move(projected.preparations.front());
  normalization_.phenotype_sum_squares=projected.response_gram.front();
  normalization_.likelihood_dimension=projected.dimension;
}
PreparedBayesianPathways::PreparedBayesianPathways(const PreparedSets& sets,
    const NamedValues& y,const BayesianPathwayOptions& o) {
  const auto n=sets.universe().ids.size(),f=sets.names().size();
  require(f>0&&f<=o.maximum_features&&n>1,"Invalid pathway dimensions");
  auto labels=sets.names();labels.push_back("(Intercept)");names(labels);
  require(y.ids==sets.universe().ids,"Response gene identities/order differ");
  require(!o.trait.empty()&&!o.response_definition.empty()&&!o.covariance_provenance.empty(),
    "Trait, response definition and covariance provenance required");
  require(std::isfinite(o.minimum_reciprocal_condition)&&o.minimum_reciprocal_condition>0&&o.minimum_reciprocal_condition<1,
    "Invalid condition threshold");
  // Two dense Gram payloads during checked binding; vectors and sparse
  // membership structures are not a process-peak estimate.
  std::size_t bytes=0;budget_add(bytes,product(2,product(f,f)),o);
  budget_add(bytes,product(4,f),o);budget_add(bytes,n,o);
  auto x=set_cross_products(sets,y,true,o.scale_features,o.scale_response,o.maximum_preparation_bytes);
  double response_max=0;for(double v:y.values)response_max=std::max(response_max,std::abs(v));
  const double relative_scale=x.response_scale/response_max;
  require(response_max>0&&x.yy*relative_scale*relative_scale>1e-24*static_cast<double>(n),
    "No usable residual response variation");
  statistics_.cross_products=gbayes::dense_cross_products(f,std::move(x.xx));
  statistics_.score=std::move(x.xy);statistics_.marker_ids=sets.names();statistics_.trait=o.trait;
  statistics_.effect_units="prepared gene-response / prepared pathway feature";
  statistics_.provenance=o.response_definition+"; "+o.covariance_provenance+
    "; direct binary-set cross-products; intercept projection of response AND pathways";
  normalization_.phenotype_sum_squares=x.yy;
  normalization_.likelihood_dimension=static_cast<double>(n-1);
  preparation_.gene_count=n;preparation_.gene_ids=sets.universe().ids;
  preparation_.nuisance_rank=1;preparation_.nuisance_names={"(Intercept)"};
  preparation_.block_count=n;preparation_.nuisance_condition=1;preparation_.worst_block_condition=1;
  preparation_.feature_scales=std::move(x.scales);preparation_.response_scale=x.response_scale;
  preparation_.residual_sum_squares=x.yy;preparation_.estimated_preparation_bytes=bytes;
  preparation_.options=o;
}
PreparedMultitraitBayesianPathways::PreparedMultitraitBayesianPathways(const CompetitiveDesign& d,
    const MultitraitGeneResponse& y,const std::vector<GeneCovarianceBlock>& blocks,
    const PathwayTraitCovariance& residual,const BayesianPathwayOptions& o) {
  const auto t=y.trait_ids.size();
  require(t>=2&&t<=64,"Multitrait pathways require 2..64 traits");
  require(o.trait.empty(),"Use response.trait_ids, leaving options.trait empty for multitrait preparation");
  require(residual.trait_ids==y.trait_ids&&residual.covariance.size()==product(t,t),"Residual covariance trait labels/order or dimensions differ");
  require(!residual.provenance.empty(),"Residual trait covariance provenance required");
  require(std::isfinite(o.minimum_reciprocal_condition)&&o.minimum_reciprocal_condition>0&&o.minimum_reciprocal_condition<1,"Invalid condition threshold");
  std::size_t covariance_bytes=0;budget_add(covariance_bytes,product(8,product(t,t)),o);
  preparation_.original_trait_covariance_condition=condition(residual.covariance,t,o.minimum_reciprocal_condition);
  auto projection_options=o;projection_options.maximum_preparation_bytes-=covariance_bytes;
  auto projected=project_pathways(d,y.gene_ids,y.values,y.trait_ids,blocks,projection_options);
  preparation_.traits=std::move(projected.preparations);
  for(auto& prepared:preparation_.traits) {
    const auto label=prepared.options.trait;prepared.options=o;prepared.options.trait=label;
    prepared.estimated_preparation_bytes+=covariance_bytes;
  }
  preparation_.response_cross_products=std::move(projected.response_gram);
  preparation_.residual_covariance=residual;
  statistics_.marker_ids=d.feature_names;statistics_.trait_ids=y.trait_ids;
  statistics_.sampling_error_covariance.resize(t*t);
  for(std::size_t j=0;j<t;++j)for(std::size_t i=0;i<t;++i)
    statistics_.sampling_error_covariance[i+j*t]=residual.covariance[i+j*t]/preparation_.traits[i].response_scale/preparation_.traits[j].response_scale;
  preparation_.prepared_trait_covariance_condition=condition(statistics_.sampling_error_covariance,t,o.minimum_reciprocal_condition);
  const auto f=d.feature_names.size();
  for(std::size_t h=0;h<t;++h) {
    gbayes::MultivariateProvider provider;provider.id=y.trait_ids[h];provider.trait=h;
    provider.statistics=std::move(projected.statistics[h]);
    provider.statistics.provenance+="; residual trait covariance: "+residual.provenance;
    provider.local_to_global.resize(f);for(std::size_t j=0;j<f;++j)provider.local_to_global[j]=j;
    provider.residual_variance=statistics_.sampling_error_covariance[h+h*t];
    // This reused field records the residual gene-response dimension, not GWAS N.
    provider.sample_size=projected.dimension;
    statistics_.providers.push_back(std::move(provider));
  }
}
MultitraitBayesianPathwayResult PreparedMultitraitBayesianPathways::fit(const gbayes::MultivariatePrior& prior,
    const gbayes::MultivariateControls& controls,const gbayes::MultivariateLearning& learning) const {
  MultitraitBayesianPathwayResult out;
  out.native=gbayes::fit_multivariate_mixture(statistics_,prior,controls,learning);out.preparation=preparation_;
  const auto t=statistics_.trait_ids.size();
  for(std::size_t j=0;j<statistics_.marker_ids.size();++j) {
    MultitraitPathwayPosterior row;row.name=statistics_.marker_ids[j];row.any_inclusion_probability=out.native.any_pip[j];
    for(std::size_t h=0;h<t;++h) {
      const auto& prepared=preparation_.traits[h];const auto scale=prepared.response_scale/prepared.feature_scales[j];
      const auto index=j*t+h;const double mean=out.native.mean[index]*scale,variance=out.native.variance[index]*scale*scale;
      require(std::isfinite(mean)&&std::isfinite(variance),"Multitrait posterior outside original-unit numerical range");
      row.mean.push_back(mean);row.variance.push_back(variance);row.inclusion_probability.push_back(out.native.pip[index]);
    }
    out.pathways.push_back(std::move(row));
  }
  return out;
}

BayesianPathwayCResult PreparedBayesianPathways::fit_bayesc(const gbayes::BayesCPrior& p,const gbayes::BayesCControls& c,gbayes::GaussianLearning l) const {
  BayesianPathwayCResult out;out.native=gbayes::fit_bayesc(statistics_,p,c,learning_for(l,normalization_));out.preparation=preparation_;out.pathways=posterior(out.native,preparation_);return out;
}
BayesianPathwayRResult PreparedBayesianPathways::fit_bayesr(const gbayes::BayesRPrior& p,const gbayes::BayesCControls& c,gbayes::GaussianLearning l) const {
  BayesianPathwayRResult out;out.native=gbayes::fit_bayesr(statistics_,p,c,learning_for(l,normalization_));out.preparation=preparation_;out.pathways=posterior(out.native,preparation_);return out;
}
}
