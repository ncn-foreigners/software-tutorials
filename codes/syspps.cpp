// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector syspps_cpp(NumericVector x, int n) {
  int N = x.size();
  IntegerVector U = sample(N, N) - 1; // sample and subtract 1 for C++ 0-based indexing
  NumericVector xx(N);
  for (int i = 0; i < N; ++i) {
    xx[i] = x[U[i]];
  }
  NumericVector z(N);
  for (int i = 0; i < N; ++i) {
    z[i] = n * sum(xx[Range(0, i)]) / sum(x);
  }
  double r = ::Rf_runif(0, 1);
  IntegerVector s(n);
  int j = 0;
  for (int i = 0; i < N; ++i) {
    if (z[i] >= r) {
      s[j] = U[i];
      r += 1;
      if (++j == n) {
        break;
      }
    }
  }
  return s.sort() + 1; // Add 1 for R's 1-based indexing
}