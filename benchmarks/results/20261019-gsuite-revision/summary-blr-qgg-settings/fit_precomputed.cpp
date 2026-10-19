#include <gbayes/bayesr.hpp>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <windows.h>
using Clock=std::chrono::steady_clock;
double elapsed(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
std::vector<double> binary(const std::string& p,std::size_t n){std::vector<double> x(n);std::ifstream f(p,std::ios::binary);f.read(reinterpret_cast<char*>(x.data()),n*sizeof(double));if(!f||f.peek()!=EOF)throw std::runtime_error("Invalid binary dimensions: "+p);return x;}
int main(int argc,char**argv){if(argc!=3)return 2;std::string in=argv[1],out=argv[2];std::ofstream(out+"/pid.txt")<<GetCurrentProcessId();
try{
 auto load=Clock::now();std::size_t n,p;double yy,ve,vb,ssb,sse;std::uint32_t seed;
 std::ifstream meta(in+"/meta.txt");meta>>n>>p>>yy>>ve>>vb>>ssb>>sse>>seed;if(!meta||p!=1324||n!=17955)throw std::runtime_error("Unexpected workload metadata");
 auto c=binary(in+"/XX.bin",p*p),xy=binary(in+"/Xy.bin",p);
 std::vector<std::string> names;std::string name;std::ifstream ids(in+"/features.txt");while(std::getline(ids,name))names.push_back(name);if(names.size()!=p)throw std::runtime_error("Feature labels differ");
 double loading=elapsed(load);auto binding=Clock::now();
 gbayes::FixedGaussianStatistics data;data.cross_products=gbayes::dense_cross_products(p,std::move(c));data.score=std::move(xy);data.marker_ids=names;data.trait="T2D gene association score";data.effect_units="Standardized gene score per standardized Reactome feature";data.provenance="Frozen installed database; sparse centered/scaled XX/Xy checked against installed qgg::computeStat";
 gbayes::BayesRPrior prior;prior.residual_variance=ve;prior.slab_variance=vb;prior.weights={.95,.02,.02,.01};prior.variance_multipliers={0,.01,.1,1};
 gbayes::BayesCControls control;control.burnin=1000;control.sampling_sweeps=5000;control.thinning=1;control.chain_seeds={seed};control.threads=1;control.retain_diagnostics=true;
 gbayes::GaussianLearning learn;learn.slab_variance=true;learn.residual_variance=true;learn.mixture_weights=true;learn.slab_prior={4,ssb};learn.residual_prior={4,sse};learn.weight_prior={1,1,1,1};learn.normalization={yy,double(n),gbayes::GaussianSpace::full};learn.residual_policy=gbayes::ResidualPolicy::exact_in_sample;
 double bind_seconds=elapsed(binding);auto start=Clock::now();
 auto fit=gbayes::fit_bayesr(data,prior,control,learn);double seconds=elapsed(start);
 if(fit.chains.size()!=1||fit.chains[0].posterior.retained_draws!=5000||fit.mean.size()!=p)throw std::runtime_error("Unexpected fit dimensions");
 auto serialization=Clock::now();std::ofstream coefficients(out+"/coefficients.csv");coefficients<<std::setprecision(17)<<"feature,mean,variance,pip\n";
 for(std::size_t j=0;j<p;++j){if(!std::isfinite(fit.mean[j])||!std::isfinite(fit.pip[j])||fit.pip[j]<0||fit.pip[j]>1)throw std::runtime_error("Nonfinite posterior");coefficients<<names[j]<<','<<fit.mean[j]<<','<<fit.variance[j]<<','<<fit.pip[j]<<'\n';}coefficients.close();
 std::ofstream trace(out+"/parameter-trace.csv");trace<<std::setprecision(17)<<"retained_draw,vb,ve,null,small,medium,large\n";const auto& t=fit.chains[0].posterior.trace.parameters;for(std::size_t i=0;i<t.size();++i){trace<<i+1<<','<<t[i].slab_variance<<','<<t[i].residual_variance;for(double w:t[i].weights)trace<<','<<w;trace<<'\n';}trace.close();
 double save_seconds=elapsed(serialization);std::ofstream result(out+"/timing.json");result<<std::setprecision(17)<<"{\n\"status\":\"completed\",\n\"statistics_loading_seconds\":"<<loading<<",\n\"statistics_binding_seconds\":"<<bind_seconds<<",\n\"sampling_seconds\":"<<seconds<<",\n\"serialization_seconds\":"<<save_seconds<<",\n\"genes\":"<<n<<",\n\"features\":"<<p<<",\n\"chains\":1,\n\"burnin\":1000,\n\"sampling_sweeps\":5000,\n\"seed\":"<<seed<<",\n\"threads_used\":"<<fit.threads_used<<",\n\"retained_draws\":"<<fit.chains[0].posterior.retained_draws<<",\n\"maximum_residual_drift\":"<<fit.chains[0].posterior.maximum_residual_drift<<"\n}\n";
 std::cout<<"Completed native fit from precomputed cross-products: "<<seconds<<" seconds\n";return 0;
}catch(const std::exception&e){std::ofstream(out+"/failure.txt")<<e.what();std::cerr<<e.what()<<'\n';return 1;}}
