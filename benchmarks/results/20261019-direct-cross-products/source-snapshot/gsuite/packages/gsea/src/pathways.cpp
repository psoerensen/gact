#include "gsea_bayes_controls.h"
#include <gbayes/diagnostics.hpp>
#include <gsea/preranked.hpp>
#include <gsea/multitrait_pathways.hpp>
#include <memory>
#include <chrono>

namespace {
using namespace Rcpp;
using Strings=std::vector<std::string>;
using Values=std::vector<double>;
gsea::Adjustment adjustment(const List& c) {
  const auto s=as<std::string>(c["adjustment"]);
  return s=="BH"?gsea::Adjustment::benjamini_hochberg:s=="BY"?gsea::Adjustment::benjamini_yekutieli:gsea::Adjustment::none;
}
NumericVector labelled(const Values& v,const Strings& names) {
  NumericVector x=wrap(v); if(!v.empty()) x.names()=wrap(names); return x;
}
NumericMatrix square(const Values& v,const Strings& names) {
  NumericMatrix x(names.size(),names.size());std::copy(v.begin(),v.end(),x.begin());
  x.attr("dimnames")=List::create(names,names);return x;
}
NumericMatrix rows(const Values& v,const Strings& row,const Strings& col) {
  if(v.empty()) return NumericMatrix(0,0);
  NumericMatrix x(row.size(),col.size());
  for(size_t i=0;i<row.size();++i)for(size_t j=0;j<col.size();++j)x(i,j)=v[i*col.size()+j];
  x.attr("dimnames")=List::create(row,col);return x;
}
List preparation(const gsea::BayesianPathwayPreparation& p,const Strings& features) {
  return List::create(_["gene_count"]=double(p.gene_count),_["nuisance_rank"]=double(p.nuisance_rank),
    _["block_count"]=double(p.block_count),_["nuisance_names"]=p.nuisance_names,
    _["feature_scales"]=labelled(p.feature_scales,features),_["response_scale"]=p.response_scale,
    _["residual_sum_squares"]=p.residual_sum_squares,_["nuisance_condition"]=p.nuisance_condition,
    _["worst_block_condition"]=p.worst_block_condition,_["estimated_preparation_bytes"]=double(p.estimated_preparation_bytes));
}
List parameters(const gbayes::GaussianParameters& p) {
  return List::create(_["residual_variance"]=p.residual_variance,_["effect_variance"]=p.slab_variance,_["weights"]=p.weights);
}
const gbayes::BayesCChain& chain(const gbayes::BayesCChain& c){return c;}
const gbayes::BayesCChain& chain(const gbayes::BayesRChain& c){return c.posterior;}
template<class Fit> List scalar(const Fit& f,const std::string& trait,double threshold) {
  const auto& n=f.native;List chains(n.chains.size());
  for(size_t k=0;k<n.chains.size();++k) {
    const auto& c=chain(n.chains[k]);List trace(c.trace.parameters.size());
    for(size_t j=0;j<c.trace.parameters.size();++j)trace[j]=parameters(c.trace.parameters[j]);
    chains[k]=List::create(_["seed"]=double(c.seed),_["retained_draws"]=double(c.retained_draws),
      _["mean"]=labelled(c.mean,n.marker_ids),_["pip"]=labelled(c.pip,n.marker_ids),
      _["parameter_mean"]=parameters(c.parameter_mean),_["trace"]=trace,
      _["maximum_residual_drift"]=c.maximum_residual_drift);
  }
  const auto a=gbayes::summarize_chains(n,threshold);Values ranges;
  for(const auto& p:a.pip)ranges.push_back(a.multiple_chains?p.maximum-p.minimum:NA_REAL);
  Strings flagged;for(auto i:a.flagged_markers)flagged.push_back(n.marker_ids[i]);
  Strings name;Values mean,var,pip;
  for(const auto& p:f.pathways){name.push_back(p.name);mean.push_back(p.mean);var.push_back(p.variance);pip.push_back(p.inclusion_probability);}
  return List::create(_["pathways"]=DataFrame::create(_["trait"]=CharacterVector(name.size(),trait),_["pathway"]=name,_["mean"]=mean,_["variance"]=var,_["pip"]=pip),
    _["parameters"]=parameters(n.parameter_mean),_["chains"]=chains,_["preparation"]=preparation(f.preparation,n.marker_ids),
    _["diagnostics"]=List::create(_["threads_used"]=n.threads_used,_["chain_pip_range"]=labelled(ranges,n.marker_ids),_["flagged_pathways"]=flagged),
    _["details"]=List::create(_["prepared_mean"]=labelled(n.mean,n.marker_ids),_["prepared_variance"]=labelled(n.variance,n.marker_ids)));
}
List multi_parameters(const gbayes::MultivariateParameters& p,const Strings& traits,const Strings& patterns,const Strings& components) {
  return List::create(_["effect_covariance"]=square(p.marker_covariance,traits),_["pattern_weights"]=labelled(p.pattern_weights,patterns),_["component_weights"]=labelled(p.component_weights,components));
}
List multitrait(const gsea::MultitraitBayesianPathwayResult& f,const List& prior,double threshold) {
  const auto& n=f.native;IntegerMatrix pm=prior["patterns"];
  const auto patterns=as<Strings>(rownames(pm));Strings components;
  for(size_t k=0;k<n.prior.component_scales.size();++k)components.push_back("component_"+std::to_string(k+1));
  List chains(n.chains.size()),prep(n.trait_ids.size());prep.attr("names")=wrap(n.trait_ids);
  for(size_t k=0;k<n.chains.size();++k) {
    const auto& c=n.chains[k];List trace(c.trace.size());
    for(size_t j=0;j<c.trace.size();++j)trace[j]=multi_parameters(c.trace[j],n.trait_ids,patterns,components);
    chains[k]=List::create(_["seed"]=double(c.seed),_["retained_draws"]=double(c.retained_draws),
      _["mean"]=rows(c.mean,n.marker_ids,n.trait_ids),_["pip"]=rows(c.pip,n.marker_ids,n.trait_ids),
      _["any_pip"]=labelled(c.any_pip,n.marker_ids),_["parameter_mean"]=multi_parameters(c.parameter_mean,n.trait_ids,patterns,components),
      _["trace"]=trace,_["maximum_residual_drift"]=c.maximum_working_residual_drift);
  }
  for(size_t t=0;t<n.trait_ids.size();++t)prep[t]=preparation(f.preparation.traits[t],n.marker_ids);
  const auto a=gbayes::summarize_chains(n,threshold);Values ranges;
  for(const auto& p:a.pip)ranges.push_back(a.multiple_chains?p.maximum-p.minimum:NA_REAL);
  Strings flagged_pathways,flagged_traits;for(auto i:a.flagged_marker_traits){flagged_pathways.push_back(n.marker_ids[i/n.trait_ids.size()]);flagged_traits.push_back(n.trait_ids[i%n.trait_ids.size()]);}
  Strings name,trait;Values mean,var,pip;
  for(size_t t=0;t<n.trait_ids.size();++t)for(const auto& p:f.pathways){name.push_back(p.name);trait.push_back(n.trait_ids[t]);mean.push_back(p.mean[t]);var.push_back(p.variance[t]);pip.push_back(p.inclusion_probability[t]);}
  return List::create(_["pathways"]=DataFrame::create(_["trait"]=trait,_["pathway"]=name,_["mean"]=mean,_["variance"]=var,_["pip"]=pip),
    _["joint"]=DataFrame::create(_["pathway"]=n.marker_ids,_["any_pip"]=n.any_pip),
    _["parameters"]=multi_parameters(n.parameter_mean,n.trait_ids,patterns,components),_["chains"]=chains,_["preparation"]=prep,
    _["diagnostics"]=List::create(_["threads_used"]=n.threads_used,_["chain_pip_range"]=rows(ranges,n.marker_ids,n.trait_ids),
      _["flagged_pathway_traits"]=DataFrame::create(_["pathway"]=flagged_pathways,_["trait"]=flagged_traits),
      _["estimated_native_memory_bytes"]=double(n.estimated_memory_bytes),_["original_trait_covariance_condition"]=f.preparation.original_trait_covariance_condition,
      _["prepared_trait_covariance_condition"]=f.preparation.prepared_trait_covariance_condition),
    _["details"]=List::create(_["prepared_mean"]=rows(n.mean,n.marker_ids,n.trait_ids),_["prepared_variance"]=rows(n.variance,n.marker_ids,n.trait_ids),
      _["pattern_probability"]=rows(n.pattern_probability,n.marker_ids,patterns),_["component_probability"]=rows(n.component_probability,n.marker_ids,components),
      _["expected_trait_sharing"]=square(n.expected_trait_sharing,n.trait_ids),_["latent_correlation"]=square(n.latent_correlation,n.trait_ids),
      _["response_cross_products"]=square(f.preparation.response_cross_products,n.trait_ids),_["prepared_residual_covariance"]=square(n.sampling_error_covariance,n.trait_ids)));
}
List competitive_result(const std::vector<gsea::CompetitiveResult>& fits,const std::string& trait) {
  Strings names;Values estimate,se,t,p,two,adjusted;List details(fits.size()),diagnostics(fits.size());
  for(size_t k=0;k<fits.size();++k) {
    const auto& f=fits[k];Strings all;Values beta;
    for(const auto& c:f.coefficients){all.push_back(c.name);beta.push_back(c.estimate);if(c.feature){names.push_back(c.name);estimate.push_back(c.estimate);se.push_back(c.standard_error);t.push_back(c.t);p.push_back(c.p_enrichment);two.push_back(c.p_two_sided);adjusted.push_back(c.adjusted_p.value_or(NA_REAL));}}
    details[k]=labelled(beta,all);
    diagnostics[k]=List::create(_["residual_df"]=double(f.residual_df),_["residual_scale"]=f.residual_scale,_["information_condition"]=f.information_condition,_["worst_block_condition"]=f.worst_block_condition,
      _["response_definition"]=f.response_definition,_["covariance_provenance"]=f.covariance_provenance);
  }
  return List::create(_["pathways"]=DataFrame::create(_["trait"]=CharacterVector(names.size(),trait),_["pathway"]=names,_["estimate"]=estimate,_["standard_error"]=se,_["t"]=t,_["p_value"]=p,_["p_two_sided"]=two,_["adjusted_p"]=adjusted),_["details"]=details,_["diagnostics"]=diagnostics);
}
}

// R conversion and orchestration only. All statistical calculations use native public APIs.
// [[Rcpp::export]]
Rcpp::List cpp_gsea(Rcpp::List stat,SEXP sets,Rcpp::NumericMatrix covariates,
    std::string method,Rcpp::List sampling,Rcpp::List prior,Rcpp::List control) {
  using namespace Rcpp;
  using Clock=std::chrono::steady_clock;
  const auto preparation_start=Clock::now();
  DataFrame first=stat[0];const auto genes=as<Strings>(first["gene"]);const auto traits=as<Strings>(stat.names());
  const bool bayes=method=="bayesc"||method=="bayesr", regression=method!="ora"&&method!="preranked";
  const bool direct_sets=bayes&&!Rf_isMatrix(sets)&&stat.size()==1&&
    as<bool>(sampling["independent"])&&covariates.ncol()==0;
  std::unique_ptr<gsea::PreparedSets> prepared;
  gsea::CompetitiveDesign design;design.gene_ids=genes;
  design.covariates=as<Values>(covariates);
  if(covariates.ncol())design.covariate_names=as<Strings>(colnames(covariates));
  if(Rf_isMatrix(sets)) {NumericMatrix x(sets);design.features=as<Values>(x);design.feature_names=as<Strings>(colnames(x));}
  else {
    List input(sets);const auto names=as<Strings>(input.names());std::vector<gsea::SetDefinition> definitions;
    for(R_xlen_t k=0;k<input.size();++k)definitions.push_back({names[k],as<Strings>(input[k])});
    gsea::Universe u;u.ids=genes;u.provenance="gsuite caller-declared complete gene background";
    prepared=std::make_unique<gsea::PreparedSets>(u,definitions);
  }
  std::vector<gsea::GeneCovarianceBlock> blocks;
  if(regression&&!direct_sets) {
    if(as<bool>(sampling["independent"]))blocks=gsea::independent_gene_covariance(genes.size());
    else for(SEXP item:as<List>(sampling["blocks"])) {List b(item);gsea::GeneCovarianceBlock block;for(int i:as<IntegerVector>(b["indices"]))block.indices.push_back(i-1);block.covariance=as<Values>(b["covariance"]);blocks.push_back(std::move(block));}
  }
  gsea::BayesianPathwayOptions bo;
  if(bayes) {
#define BOPT(name,type) bo.name=as<type>(control[#name])
    BOPT(scale_response,bool);BOPT(scale_features,bool);BOPT(maximum_features,size_t);BOPT(maximum_covariates,size_t);BOPT(maximum_block_size,size_t);BOPT(maximum_preparation_bytes,size_t);BOPT(minimum_reciprocal_condition,double);
#undef BOPT
    bo.response_definition="Caller-supplied gene score";bo.covariance_provenance=as<std::string>(sampling["gene_provenance"]);
    if(prepared&&!direct_sets)design=gsea::bayesian_set_design(*prepared,design.covariate_names,design.covariates,bo);
  } else if(regression&&prepared&&method!="magma")design=gsea::set_design(*prepared,design.covariate_names,design.covariates);
  if(bayes&&stat.size()>1) {
    gsea::MultitraitGeneResponse response;response.gene_ids=genes;response.trait_ids=traits;
    for(SEXP s:stat){DataFrame d(s);const auto y=as<Values>(d["score"]);response.values.insert(response.values.end(),y.begin(),y.end());}
    gsea::PathwayTraitCovariance residual{traits,as<Values>(sampling["trait_covariance"]),as<std::string>(sampling["trait_provenance"])};
    gbayes::MultivariatePrior p;p.marker_covariance=as<Values>(prior["effect_covariance"]);IntegerMatrix patterns=prior["patterns"];
    for(int k=0;k<patterns.nrow();++k){gbayes::TraitPattern pattern;for(int t=0;t<patterns.ncol();++t)pattern.push_back(patterns(k,t));p.patterns.push_back(pattern);}
    p.pattern_weights=as<Values>(prior["pattern_weights"]);p.component_scales=as<Values>(prior["component_scales"]);p.component_weights=as<Values>(prior["component_weights"]);
    gbayes::MultivariateLearning l;l.marker_covariance=as<bool>(control["estimate_covariance"]);l.pattern_weights=as<bool>(control["estimate_pattern_weights"]);l.component_weights=as<bool>(control["estimate_component_weights"]);
    if(l.marker_covariance){List iw=prior["covariance_prior"];l.covariance_prior_df=as<double>(iw["df"]);l.covariance_prior_scale=as<Values>(iw["scale"]);}
    if(l.pattern_weights)l.pattern_prior_counts=as<Values>(prior["pattern_prior"]);
    if(l.component_weights)l.component_prior_counts=as<Values>(prior["component_prior"]);
    gbayes::MultivariateControls c;c.chains=gsea_r::chain_controls(control);c.memory_limit_bytes=as<size_t>(control["memory_limit_bytes"]);c.retain_pattern_probabilities=as<bool>(control["retain_pattern_probabilities"]);
    const gsea::PreparedMultitraitBayesianPathways data(design,response,blocks,residual,bo);
    const auto ready=Clock::now();const auto fitted=data.fit(p,c,l);const auto sampled=Clock::now();
    auto out=multitrait(fitted,prior,as<double>(control["pip_range_threshold"]));
    out["timing"]=List::create(_["preparation_seconds"]=std::chrono::duration<double>(ready-preparation_start).count(),
      _["blr_fit_seconds"]=std::chrono::duration<double>(sampled-ready).count(),
      _["posterior_summary_seconds"]=std::chrono::duration<double>(Clock::now()-sampled).count(),
      _["preparation_route"]="projected_design");
    return List::create(_["fits"]=List::create(_["joint"]=out),_["native_version"]=gsea::version_string(),_["joint_traits"]=true);
  }
  List fits(stat.size());fits.attr("names")=wrap(traits);
  for(R_xlen_t k=0;k<stat.size();++k) {
    checkUserInterrupt();DataFrame s=stat[k];List out;
    if(method=="ora") {
      Strings selected;LogicalVector choose=s["selected"];for(size_t i=0;i<genes.size();++i)if(choose[i])selected.push_back(genes[i]);
      const auto f=gsea::overrepresentation(*prepared,selected,adjustment(control));Strings name;Values size,count,expected,fold,p,logp,adj;
      for(const auto& r:f.rows){name.push_back(r.name);size.push_back(r.size);count.push_back(r.selected_in_set);expected.push_back(r.expected);fold.push_back(r.fold_enrichment.value_or(NA_REAL));p.push_back(r.probability.p_value);logp.push_back(r.probability.log_p_value);adj.push_back(r.adjusted_p.value_or(NA_REAL));}
      out=List::create(_["pathways"]=DataFrame::create(_["trait"]=CharacterVector(name.size(),traits[k]),_["pathway"]=name,_["size"]=size,_["selected"]=count,_["expected"]=expected,_["fold_enrichment"]=fold,_["p_value"]=p,_["log_p_value"]=logp,_["adjusted_p"]=adj),_["diagnostics"]=List::create(_["selected_count"]=double(f.selected_count)));
    } else if(method=="preranked") {
      gsea::PrerankedOptions c;
#define POPT(name,type) c.name=as<type>(control[#name])
      POPT(exponent,double);POPT(minimum_size,size_t);POPT(maximum_size,size_t);POPT(replicates,uint64_t);POPT(seed,uint64_t);POPT(threads,unsigned);POPT(maximum_scratch_bytes,size_t);
#undef POPT
      auto dir=as<std::string>(control["direction"]);c.direction=dir=="positive"?gsea::EnrichmentDirection::positive:dir=="negative"?gsea::EnrichmentDirection::negative:gsea::EnrichmentDirection::standard;c.adjustment=adjustment(control);
      const auto f=gsea::preranked(*prepared,{genes,as<Values>(s["score"])},c);Strings name,diagnostic;Values size,es,nes,p,adj,nulls;LogicalVector tested,available;List edges;
      for(const auto& r:f.rows){name.push_back(r.name);diagnostic.push_back(r.diagnostic);size.push_back(r.size);es.push_back(r.es);nes.push_back(r.nes);p.push_back(r.p_value);adj.push_back(r.adjusted_p.value_or(NA_REAL));nulls.push_back(r.same_sign_nulls);tested.push_back(r.tested);available.push_back(r.inference_available);edges.push_back(wrap(r.leading_edge));}
      edges.attr("names")=wrap(name);
      out=List::create(_["pathways"]=DataFrame::create(_["trait"]=CharacterVector(name.size(),traits[k]),_["pathway"]=name,_["size"]=size,_["es"]=es,_["nes"]=nes,_["p_value"]=p,_["adjusted_p"]=adj,_["tested"]=tested,_["inference_available"]=available,_["same_sign_nulls"]=nulls,_["diagnostic"]=diagnostic),
        _["details"]=List::create(_["leading_edge"]=edges,_["ranking"]=labelled(f.ranking.values,f.ranking.ids)),_["diagnostics"]=List::create(_["threads_used"]=f.workers_used,_["scratch_bytes"]=double(f.scratch_bytes),_["tied_values"]=double(f.tied_values),_["tested_family_size"]=double(f.tested_family_size),_["null_model"]=f.null_model));
    } else if(method=="magma") {
      gsea::SummaryEnrichmentOptions c;
#define MOPT(name,type) c.name=as<type>(control[#name])
      MOPT(correct_size,bool);MOPT(correct_density,bool);MOPT(correct_sample_size,bool);MOPT(correct_mac,bool);MOPT(truncate,bool);MOPT(residualize_internal,bool);MOPT(lower_z_bound,double);MOPT(upper_z_sd,double);MOPT(covariate_sd_bound,double);MOPT(correlation_power,double);
#undef MOPT
      auto p=as<Values>(s["p_value"]),m=as<Values>(s["marker_count"]),e=as<Values>(s["effective_parameters"]),n=as<Values>(s["sample_size"]),a=as<Values>(s["mean_minor_allele_count"]);std::vector<gsea::GeneSummary> summaries;
      for(size_t i=0;i<genes.size();++i)summaries.push_back({genes[i],p[i],m[i],e[i],n[i],a[i]});
      const auto prep=gsea::prepare_summary_enrichment(*prepared,summaries,blocks,c);
      out=competitive_result(gsea::fit_summary_enrichment(prep,as<bool>(control["joint"]),adjustment(control)),traits[k]);
      out["preparation"]=List::create(_["response"]=labelled(prep.response.values,genes),_["transformed_response"]=labelled(prep.transformed_response.values,genes),_["constant_covariates"]=prep.constant_covariates,_["clipped_responses"]=double(prep.clipped_responses),_["clipped_covariates"]=double(prep.clipped_covariates));
    } else if(method=="competitive") {
      gsea::CompetitiveOptions c;c.maximum_block_size=as<size_t>(control["maximum_block_size"]);c.maximum_coefficients=as<size_t>(control["maximum_coefficients"]);c.minimum_reciprocal_condition=as<double>(control["minimum_reciprocal_condition"]);c.response_definition="Caller-supplied gene score";c.covariance_provenance=as<std::string>(sampling["gene_provenance"]);
      gsea::NamedValues response{genes,as<Values>(s["score"])};std::vector<gsea::CompetitiveResult> results;
      if(as<bool>(control["joint"])) results.push_back(gsea::PreparedCompetitiveModel(design,blocks,c).fit(response,adjustment(control)));
      else {
        Values ps;for(size_t j=0;j<design.feature_names.size();++j){auto d=design;d.feature_names={design.feature_names[j]};d.features.assign(design.features.begin()+j*genes.size(),design.features.begin()+(j+1)*genes.size());results.push_back(gsea::PreparedCompetitiveModel(d,blocks,c).fit(response));ps.push_back(results.back().coefficients.back().p_enrichment);}
        const auto adj=gsea::adjust_p_values(ps,adjustment(control));if(adjustment(control)!=gsea::Adjustment::none)for(size_t j=0;j<results.size();++j)results[j].coefficients.back().adjusted_p=adj[j];
      }
      out=competitive_result(results,traits[k]);
    } else {
      bo.trait=traits[k];
      const gsea::NamedValues response{genes,as<Values>(s["score"])};
      const auto data=direct_sets ?
        std::make_unique<gsea::PreparedBayesianPathways>(*prepared,response,bo) :
        std::make_unique<gsea::PreparedBayesianPathways>(design,response,blocks,bo);
      gbayes::GaussianLearning l;l.slab_variance=as<bool>(control["estimate_effect_variance"]);l.residual_variance=as<bool>(control["estimate_residual_variance"]);l.mixture_weights=as<bool>(control["estimate_weights"]);
      if(l.slab_variance){List p=prior["effect_variance_prior"];l.slab_prior={as<double>(p["df"]),as<double>(p["scale"])};}
      if(l.residual_variance){List p=prior["residual_variance_prior"];l.residual_prior={as<double>(p["df"]),as<double>(p["scale"])};}
      if(l.mixture_weights)l.weight_prior=as<Values>(prior["weight_prior"]);
      const auto c=gsea_r::chain_controls(control);const auto threshold=as<double>(control["pip_range_threshold"]);
      const auto ready=Clock::now();auto sampled=ready;
      if(method=="bayesc") {gbayes::BayesCPrior p;p.residual_variance=as<double>(prior["residual_variance"]);p.slab_variance=as<double>(prior["effect_variance"]);p.inclusion_probability=as<double>(prior["inclusion_probability"]);const auto fitted=data->fit_bayesc(p,c,l);sampled=Clock::now();out=scalar(fitted,traits[k],threshold);}
      else {gbayes::BayesRPrior p;p.residual_variance=as<double>(prior["residual_variance"]);p.slab_variance=as<double>(prior["effect_variance"]);p.weights=as<Values>(prior["weights"]);p.variance_multipliers=as<Values>(prior["variance_multipliers"]);const auto fitted=data->fit_bayesr(p,c,l);sampled=Clock::now();out=scalar(fitted,traits[k],threshold);}
      out["timing"]=List::create(_["preparation_seconds"]=std::chrono::duration<double>(ready-preparation_start).count(),
        _["blr_fit_seconds"]=std::chrono::duration<double>(sampled-ready).count(),
        _["posterior_summary_seconds"]=std::chrono::duration<double>(Clock::now()-sampled).count(),
        _["preparation_route"]=direct_sets?"direct_set_cross_products":"projected_design");
    }
    fits[k]=out;
  }
  return List::create(_["fits"]=fits,_["native_version"]=gsea::version_string(),_["joint_traits"]=false);
}


