#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector twoSidedP(NumericVector p) {
  // calculate the two sided P values (P(<=the observed probability) 
  // based on the sorted probability values p in the ascending order
  
  int n = p.size();
  int upperBoundPointer = 0;
  double accumulatedSum = 0;
  // double relErr = 1 + 1e-7;
  NumericVector twosided(n);
  for (int i = 0; i < n; i++) {
    // double upperBound = p[i] * relErr;
    double upperBound = p[i];
    for (int j = upperBoundPointer; j < n; j++) {
      if (p[j] <= upperBound) {
        accumulatedSum += p[j];
      } else {
        // update upper bound pointer
        upperBoundPointer = j;
        break;
      }
    }
    // assign the sum of all probabilities <= p[i] * relErr
    twosided[i] = accumulatedSum;
  }
  
  return(twosided);
}