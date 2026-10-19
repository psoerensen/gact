#include "gsea/set_cross_products.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_map>
namespace gsea {
SetCrossProducts set_cross_products(const PreparedSets& p,const NamedValues& y,
    bool center,bool scale_features,bool scale_response,std::size_t budget) {
  const auto n=p.universe().ids.size(),m=p.names().size();
  if(n<2||m==0) throw std::invalid_argument("Need at least two observations and one set");
  if(m>std::numeric_limits<std::size_t>::max()/m || m*m>budget/sizeof(double))
    throw std::invalid_argument("Set cross-product matrix exceeds maximum_matrix_bytes");
  if(y.ids.size()!=n||y.values.size()!=n)
    throw std::invalid_argument("Response must cover the full universe exactly once");
  std::unordered_map<std::string,double> by_id;
  for(std::size_t i=0;i<n;++i)
    if(!std::isfinite(y.values[i])||!by_id.emplace(y.ids[i],y.values[i]).second)
      throw std::invalid_argument("Response IDs must be unique and values finite");
  std::vector<double> response(n);double response_max=0;
  for(std::size_t i=0;i<n;++i) {
    auto at=by_id.find(p.universe().ids[i]);
    if(at==by_id.end()) throw std::invalid_argument("Response IDs differ from universe");
    response[i]=at->second;response_max=std::max(response_max,std::abs(response[i]));
  }
  // Protect response products before centering, as the general GLS path does.
  // An online mean also makes exactly constant responses exactly constant.
  if(response_max==0) response_max=1;
  double mean=0;
  for(std::size_t i=0;i<n;++i) {
    response[i]/=response_max;
    mean+=(response[i]-mean)/static_cast<double>(i+1);
  }
  SetCrossProducts out;out.n=n;out.feature_names=p.names();
  const double response_center=center?mean:0;
  out.response_center=response_center*response_max;
  double response_divisor=1;
  if(scale_response) {
    // Match scale(y, center=center): RMS when not centering, sample SD otherwise.
    double ss=0;for(double v:response) ss+=(v-response_center)*(v-response_center);
    response_divisor=std::sqrt(ss/static_cast<double>(n-1));
    out.response_scale=response_divisor*response_max;
    if(!(out.response_scale>0)||!std::isfinite(out.response_scale))
      throw std::invalid_argument("Cannot scale a zero or nonfinite response");
  }
  double sum_y=0;
  for(auto& v:response) {
    v=scale_response?(v-response_center)/response_divisor:(v-response_center)*response_max;
    sum_y+=v;out.yy+=v*v;
  }
  out.xx.assign(m*m,0);out.xy.assign(m,0);out.centers.resize(m);out.scales.assign(m,1);
  // Inverted membership lists accumulate only co-occurring set pairs.
  // Cost before dense centering is sum_i degree(i)^2, not a scan of every
  // pair of sets. No n-by-m feature matrix is built.
  std::vector<std::vector<std::size_t>> memberships(n);
  std::vector<double> sizes(m);
  const double dn=static_cast<double>(n);
  for(std::size_t j=0;j<m;++j) {
    const double k=static_cast<double>(p.offsets()[j+1]-p.offsets()[j]);sizes[j]=k;
    if(k==0||k==dn) throw std::invalid_argument("Constant set feature: "+p.names()[j]);
    out.centers[j]=center?k/dn:0;
    if(scale_features) out.scales[j]=std::sqrt((center?k-k*k/dn:k)/static_cast<double>(n-1));
    for(auto offset=p.offsets()[j];offset<p.offsets()[j+1];++offset) {
      const auto i=p.members()[offset];out.xy[j]+=response[i];memberships[i].push_back(j);
    }
    out.xy[j]=(out.xy[j]-out.centers[j]*sum_y)/out.scales[j];
  }
  for(const auto& a:memberships)
    for(std::size_t j=0;j<a.size();++j) for(std::size_t k=0;k<=j;++k)
      ++out.xx[a[j]*m+a[k]];
  for(std::size_t j=0;j<m;++j) for(std::size_t k=0;k<=j;++k) {
    const double v=(out.xx[j*m+k]-(center?sizes[j]*sizes[k]/dn:0))/(out.scales[j]*out.scales[k]);
    out.xx[j*m+k]=out.xx[k*m+j]=v;
  }
  if(!std::isfinite(out.yy)||!std::isfinite(sum_y)) throw std::invalid_argument("Response summaries overflow");
  for(double v:out.xy) if(!std::isfinite(v)) throw std::invalid_argument("Response cross-products overflow");
  return out;
}
}
