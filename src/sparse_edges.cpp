#include <Rcpp.h>
#include <algorithm>
#include <numeric>
using namespace Rcpp;

// Enumerate only edges with cost < 1; sorted positive masses required.
// [[Rcpp::export]]
List sparse_edges(NumericVector a, NumericVector b, double width) {
  int n=a.size(), m=b.size();
  std::vector<int> parent(n+m); std::iota(parent.begin(),parent.end(),0);
  auto root=[&](int x) { while(parent[x]!=x) {parent[x]=parent[parent[x]]; x=parent[x];} return x; };
  std::vector<int> ii,jj; std::vector<double> cc;
  double t=width/2e6;
  for(int i=0;i<n;++i) {
    // Expand bounds slightly; exact original cost below decides membership.
    double lo=t<1 ? a[i]*(1-t)/(1+t) : 0;
    double hi=t<1 ? a[i]*(1+t)/(1-t) : R_PosInf;
    auto begin=std::lower_bound(b.begin(),b.end(),lo*(1-1e-12));
    for(auto it=begin;it!=b.end() && *it<=hi*(1+1e-12);++it) {
      int j=it-b.begin();
      double c=std::abs(a[i]-b[j])/((a[i]+b[j])/2)*1e6/width;
      if(c<1) { ii.push_back(i+1);jj.push_back(j+1);cc.push_back(c); parent[root(i)]=root(n+j); }
    }
  }
  std::vector<int> component(ii.size());
  for(size_t k=0;k<ii.size();++k) component[k]=root(ii[k]-1)+1;
  return List::create(_["i"]=ii,_["j"]=jj,_["cost"]=cc,_["component"]=component);
}
