#include <Rcpp.h>
using namespace Rcpp;
// Independent sum of Euclidean edge lengths over the saved fixed routes.
// [[Rcpp::export]]
NumericVector pilot_path_lengths(NumericMatrix z, IntegerVector offsets,
                                 IntegerVector u, IntegerVector v) {
  NumericVector result(offsets.size()-1);
  for(int p=0;p<result.size();++p) {
    double length=0;
    for(int e=offsets[p];e<offsets[p+1];++e) {
      double squared=0;
      for(int d=0;d<z.ncol();++d) {
        double diff=z(u[e],d)-z(v[e],d); squared+=diff*diff;
      }
      length+=std::sqrt(squared);
    }
    result[p]=length;
  }
  return result;
}
