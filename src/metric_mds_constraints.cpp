#include <Rcpp.h>
#include <vector>
#include <numeric>
// Linear-time graph metadata for already validated, unique, lexicographic pairs.
// No path distances are computed: supplied targets remain fitting constraints.
// [[Rcpp::export]]
Rcpp::List grip_constraint_graph_cpp(int n, Rcpp::IntegerMatrix pairs, Rcpp::NumericVector targets) {
  if(n<2 || pairs.ncol()!=2 || pairs.nrow()!=targets.size())Rcpp::stop("Invalid constraint graph");
  std::vector<int> degree(n,0),offset(n,0),parent(n),rank(n,0);
  std::iota(parent.begin(),parent.end(),0);
  auto find=[&](int a){while(parent[a]!=a){parent[a]=parent[parent[a]];a=parent[a];}return a;};
  for(int q=0;q<pairs.nrow();++q){
    int a=pairs(q,0)-1,b=pairs(q,1)-1;
    if(a<0||b<0||a>=n||b>=n||a==b||!R_FINITE(targets[q])||targets[q]<=0)Rcpp::stop("Invalid constraint edge");
    ++degree[a];++degree[b];int ra=find(a),rb=find(b);
    if(ra!=rb){if(rank[ra]<rank[rb])std::swap(ra,rb);parent[rb]=ra;if(rank[ra]==rank[rb])++rank[ra];}
  }
  Rcpp::List adj(n),weight(n);
  for(int i=0;i<n;++i){adj[i]=Rcpp::IntegerVector(degree[i]);weight[i]=Rcpp::NumericVector(degree[i]);}
  for(int q=0;q<pairs.nrow();++q){
    if((q&1048575)==0)Rcpp::checkUserInterrupt();
    int a=pairs(q,0)-1,b=pairs(q,1)-1;
    Rcpp::IntegerVector aa=adj[a],ab=adj[b];Rcpp::NumericVector wa=weight[a],wb=weight[b];
    aa[offset[a]]=b+1;wa[offset[a]++]=targets[q];ab[offset[b]]=a+1;wb[offset[b]++]=targets[q];
  }
  Rcpp::IntegerVector component(n);std::vector<int> labels(n,0);int count=0;
  for(int i=0;i<n;++i){int r=find(i);if(!labels[r])labels[r]=++count;component[i]=labels[r];}
  return Rcpp::List::create(Rcpp::_["adj_list"]=adj,Rcpp::_["weight_list"]=weight,
    Rcpp::_["component_id"]=component,Rcpp::_["n_components"]=count);
}
