#include <gsea/summary_enrichment.hpp>
#include <boost/math/distributions/normal.hpp>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <map>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <stdexcept>
using Clock=std::chrono::steady_clock;
std::vector<std::string> split(const std::string& s){std::istringstream in(s);std::vector<std::string> v;std::string t;while(in>>t)v.push_back(t);return v;}
int main(int argc,char** argv){double loading=0,preparation=0,fitting=0;std::string phase="input_loading";auto phase_start=Clock::now();try{
 if(argc!=5)throw std::runtime_error("usage: classical_driver genes.raw genes.out.txt sets.tsv output_prefix");
 std::ifstream pin(argv[2]);std::string line;std::unordered_map<std::string,double> p;std::size_t pi=0;bool header=false;
 while(std::getline(pin,line)){if(line.empty()||line[0]=='#')continue;auto v=split(line);if(!header){auto it=std::find(v.begin(),v.end(),"P");if(it==v.end())throw std::runtime_error("No P header");pi=it-v.begin();header=true;}else p[v[0]]=std::stod(v.at(pi));}
 std::ifstream raw(argv[1]);std::vector<gsea::GeneSummary> rows;std::vector<std::vector<double>> suffix;gsea::Universe universe;std::vector<std::size_t> parent;
 while(std::getline(raw,line)){if(line.empty()||line[0]=='#')continue;auto v=split(line);if(v.size()<9)throw std::runtime_error("Invalid MAGMA raw row");auto i=rows.size();rows.push_back({v[0],boost::math::cdf(boost::math::complement(boost::math::normal(),std::stod(v[8]))),std::stod(v[4]),std::stod(v[5]),std::stod(v[6]),std::stod(v[7])});universe.ids.push_back(v[0]);parent.push_back(i);std::vector<double> c;for(std::size_t k=9;k<v.size();++k)c.push_back(std::stod(v[k]));if(c.size()>i)throw std::runtime_error("Invalid raw correlation suffix");suffix.push_back(c);}
 auto root=[&](std::size_t i){while(parent[i]!=i)i=parent[i];return i;};
 for(std::size_t i=0;i<rows.size();++i)for(std::size_t j=0;j<suffix[i].size();++j)if(suffix[i][j]!=0)parent[root(i)]=root(i-suffix[i].size()+j);
 std::map<std::size_t,std::vector<std::size_t>> groups;for(std::size_t i=0;i<rows.size();++i)groups[root(i)].push_back(i);
 std::vector<gsea::GeneCovarianceBlock> blocks;std::vector<std::size_t> bi(rows.size()),local(rows.size());
 for(auto& pair:groups){auto ids=pair.second;if(ids.size()>2048)throw std::runtime_error("Gene-correlation component exceeds installed 2048-gene bound; no silent split");gsea::GeneCovarianceBlock b{ids,std::vector<double>(ids.size()*ids.size(),0)};for(std::size_t j=0;j<ids.size();++j){b.covariance[j+j*ids.size()]=1;bi[ids[j]]=blocks.size();local[ids[j]]=j;}blocks.push_back(std::move(b));}
 for(std::size_t i=0;i<rows.size();++i)for(std::size_t j=0;j<suffix[i].size();++j){auto other=i-suffix[i].size()+j;double r=suffix[i][j];if(r==0)continue;auto& b=blocks[bi[i]];auto n=b.indices.size();b.covariance[local[i]+local[other]*n]=r;b.covariance[local[other]+local[i]*n]=r;}
 std::unordered_map<std::string,bool> known;for(auto& id:universe.ids)known[id]=true;
 std::ifstream sin(argv[3]);std::vector<gsea::SetDefinition> definitions;
 while(std::getline(sin,line)){auto v=split(line);if(v.empty())continue;gsea::SetDefinition s{v[0],{}};for(std::size_t j=1;j<v.size();++j){if(!known.count(v[j]))throw std::runtime_error("Gene universe differs from frozen sets");s.members.push_back(v[j]);}definitions.push_back(std::move(s));}
 loading=std::chrono::duration<double>(Clock::now()-phase_start).count();
 std::cout<<rows.size()<<" genes, "<<definitions.size()<<" pathways, "<<blocks.size()<<" correlation blocks loaded\n"<<std::flush;
 std::size_t rank=0,shared_applications=0,quadratic_entries=0;
 std::ofstream result(std::string(argv[4])+".csv");result<<std::setprecision(17)<<"set,beta,se,p\n";
 for(std::size_t start=0;start<definitions.size();start+=definitions.size()){auto end=std::min(start+definitions.size(),definitions.size());phase="preparation";phase_start=Clock::now();auto t=phase_start;gsea::PreparedSets sets(universe,std::vector<gsea::SetDefinition>(definitions.begin()+start,definitions.begin()+end));auto prepared=gsea::prepare_summary_enrichment(sets,rows,blocks);rank=prepared.prepared_covariance->rank();shared_applications+=prepared.marginal_cross_products->shared_block_applications();quadratic_entries+=prepared.marginal_cross_products->membership_quadratic_entries();preparation+=std::chrono::duration<double>(Clock::now()-t).count();phase="fitting";phase_start=Clock::now();t=phase_start;auto fits=gsea::fit_summary_enrichment(prepared);fitting+=std::chrono::duration<double>(Clock::now()-t).count();for(auto& fit:fits){auto c=fit.coefficients.back();result<<c.name<<','<<c.estimate<<','<<c.standard_error<<','<<c.p_enrichment<<'\n';}}
 std::ofstream timing(std::string(argv[4])+".timing.json");timing<<std::setprecision(17)<<"{\"input_loading_seconds\":"<<loading<<",\"covariance_rank\":"<<rank<<",\"shared_block_applications\":"<<shared_applications<<",\"membership_quadratic_entries\":"<<quadratic_entries<<",\"preparation_seconds\":"<<preparation<<",\"classical_fit_seconds\":"<<fitting<<",\"genes\":"<<rows.size()<<",\"pathways\":"<<definitions.size()<<",\"batch_size\":1324,\"shared_gene_input\":true}";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';if(argc==5){std::ofstream failure(std::string(argv[4])+".failure.json");failure<<std::setprecision(17)<<"{\"status\":\"failed\",\"phase\":\""<<phase<<"\",\"input_loading_seconds\":"<<loading<<",\"completed_preparation_seconds\":"<<preparation<<",\"completed_fitting_seconds\":"<<fitting<<",\"failed_phase_elapsed_seconds\":"<<std::chrono::duration<double>(Clock::now()-phase_start).count()<<"}";}return 1;}}

