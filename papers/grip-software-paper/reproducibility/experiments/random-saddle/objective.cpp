#include <Rcpp.h>
using namespace Rcpp;

// Exact squared profiled relative RMSE, without a length stabilizer.
// path[p] = sum_{e in gamma[p]} ||z_u-z_v||; scale = <path,g>/<g,g>.
// [[Rcpp::export]]
List saddle_objective(NumericVector coords, int n, IntegerMatrix edges,
                      IntegerVector offsets, IntegerVector path_edges,
                      NumericVector targets, bool edge_only = false) {
  int dim = coords.size()/n, m = edges.nrow();
  NumericVector lengths(m), edge_gradient(m), gradient(coords.size());
  for (int e=0;e<m;++e) {
    double d2=0;
    for (int d=0;d<dim;++d) {
      double v=coords[edges(e,0)-1+n*d]-coords[edges(e,1)-1+n*d]; d2+=v*v;
    }
    lengths[e]=std::sqrt(d2);
  }
  int np = targets.size();
  NumericVector observed(np);
  double gg=0, yg=0, yy=0;
  for (int p=0;p<np;++p) {
    double y=0;
    if (edge_only) y=lengths[p];
    else for(int j=offsets[p];j<offsets[p+1];++j) y+=lengths[path_edges[j]];
    observed[p]=y; gg+=targets[p]*targets[p]; yg+=y*targets[p]; yy+=y*y;
  }
  if (!(yg>1e-50) || !(gg>0)) stop("Degenerate coordinate scale");
  double scale=yg/gg, residual=0;
  for (int p=0;p<np;++p) {
    double r=observed[p]-scale*targets[p]; residual+=r*r;
    double dp=2*gg/(yg*yg)*(observed[p]-yy/yg*targets[p]);
    if (edge_only) edge_gradient[p]=dp;
    else for(int j=offsets[p];j<offsets[p+1];++j) edge_gradient[path_edges[j]]+=dp;
  }
  for (int e=0;e<m;++e) if(lengths[e]>0) {
    for(int d=0;d<dim;++d) {
      int u=edges(e,0)-1+n*d, v=edges(e,1)-1+n*d;
      double grad=edge_gradient[e]*(coords[u]-coords[v])/lengths[e];
      gradient[u]+=grad; gradient[v]-=grad;
    }
  }
  return List::create(_["value"]=residual/(scale*scale*gg),
                      _["gradient"]=gradient, _["scale"]=scale);
}
