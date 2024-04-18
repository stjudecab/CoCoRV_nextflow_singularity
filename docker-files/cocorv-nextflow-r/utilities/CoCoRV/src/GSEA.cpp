// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <algorithm>    // std::sort, std::stable_sort
#include <string>
#include <math.h>

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace std;

// sort using the decreasing order and return the index
IntegerVector sort_indexes_decreasing(NumericVector &v) {
  
  // initialize original index locations
  IntegerVector idx = seq_along(v) - 1;
  
  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  stable_sort(idx.begin(), idx.end(),
              [&v](int i1, int i2) {return v[i1] > v[i2];});
  
  return idx;
}

// [[Rcpp::export]]
double calculateESCpp(NumericVector score, LogicalVector inSet, std::string scoreType) {
  // calculate the ES of GSEA from the score
  
  // input
  // score: the score of each gene, sorted in a decreasing order
  // inSet: the indicator of being in a set
  
  // output
  // ES: the ES defined in GSEA
  
  int k = 0; 
  double NS = 0;
  int N = score.length();
  for (int i = 0; i < N; i++) {
    if (inSet[i]) {
      k++;
      NS += abs(score[i]);
    }
  }
  
  // cout<<"k="<<k<<"NS="<<NS<<"N="<<N<<endl;
  
  NumericVector ESi(N);
  double delta;
  int maxIndex = 0;
  int minIndex = 0;
  double splus = -numeric_limits<double>::max();
  double sminus = numeric_limits<double>::max();
  for (int i = 0; i < N; i++) {
    if (inSet[i]) {
      if (NS > 0) {
        delta = abs(score[i]) / NS;
      } else {
        // treat each gene equally
        delta = 1.0 / k;
      }
    } else {
      delta = - 1.0 / (N - k);
    }
    if (i == 0) {
      ESi[0] = delta;
    } else {
      ESi[i] = ESi[i - 1] + delta;
    }
    if (ESi[i] > splus) {
      splus = ESi[i];
      maxIndex = i;
    }
    if (ESi[i] < sminus) {
      sminus = ESi[i];
      minIndex = i;
    }
    // cout << "score[i]" << score[i] << "delta=" << delta <<
    // " ESi="<< ESi[i] << " splus=" << splus << " sminus=" << sminus << endl;
  }
  
  if (scoreType == "std") {
    if (abs(splus) >= abs(sminus)) {
      return splus;
    } else {
      return sminus;
    }
  } else if (scoreType == "pos") {
    return splus;
  } else if (scoreType == "neg") {
    return sminus;
  }
}

// [[Rcpp::export]]
NumericVector calculateESCppFast(NumericVector score, LogicalVector inSet, std::string scoreType, double p, bool returnLeadingEdge) {
  // calculate the ES of GSEA from the score
  
  // input
  // score: the score of each gene, sorted in a decreasing order
  // inSet: the indicator of being in a set
  // scoreType: the score type
  // p: the power parameter
  
  // output
  // ES: the ES defined in GSEA
  
  int k = 0; 
  int N = score.length();
  double NS = 0;
  
  for (int i = 0; i < N; i++) {
    if (inSet[i]) {
      k++;
    }
  }
  
  IntegerVector selected(k);
  NumericVector cumulativeScores(k); 
  NumericVector adjustedScores(k);
  
  int j = 0;
  for (int i = 0; i < N; i++) {
    if (inSet[i]) {
      double adjusted = pow(abs(score[i]), p);
      NS += adjusted;
      adjustedScores[j] = adjusted;
      selected[j] = i;
      cumulativeScores[j] = NS;
      j++;
    }
  }
  
  // cout<<"k="<<k<<"NS="<<NS<<"N="<<N<<endl;
  
  NumericVector ESiMax(k);  // right at one gene in the pathway
  NumericVector ESiMin(k);  // right before one gene in the pathway
  double delta;
  int maxIndex = 0;
  int minIndex = 0;
  double splus = -numeric_limits<double>::max();
  double sminus = numeric_limits<double>::max();
  for (int i = 0; i < k; i++) {
    double firstPart = 0;
    if (NS > 0) {
      firstPart = cumulativeScores[i] / NS;
    } else {
      firstPart = (i + 1.0) / k;
    }
    ESiMax[i] = firstPart - (selected[i] - i) * 1.0 / (N - k); 
    if (NS > 0) {
      // the score before the gene
      ESiMin[i] = ESiMax[i] - adjustedScores[i] / NS;  
    } else {
      ESiMin[i] = ESiMax[i] - 1.0 / k;
    }
    if (ESiMax[i] > splus) {
      splus = ESiMax[i];
      maxIndex = i;
    }
    if (ESiMin[i] < sminus) {
      sminus = ESiMin[i];
      minIndex = i;
    }
     // cout << "adjustedScores[i]" << adjustedScores[i] << "ESiMax[i]"
     //      << ESiMax[i] << " ESiMin="<< ESiMin[i] << " splus=" << splus
     //      << " sminus=" << sminus << "maxIndex=" << maxIndex << 
     //      "minIndex" << minIndex << endl;
  }
  
  
  
  double ES = 0;
  if (scoreType == "std") {
    if (splus >= -sminus) {
      ES = splus;
    } else {
      ES = sminus;
    }
  } else if (scoreType == "pos") {
    ES = splus;
  } else if (scoreType == "neg") {
    ES = sminus;
  }
  
  NumericVector returnValues = {ES};

  if (returnLeadingEdge) {
    if (scoreType == "std") {
      if (splus >= -sminus) {
        for (int i = 0; i <= maxIndex; i++) {
          returnValues.push_back(selected[i]);
        }
      } else {
        for (int i = k - 1; i >= minIndex; i--) {
          returnValues.push_back(selected[i]);
        }
      }
    } else if (scoreType == "pos") {
      for (int i = 0; i <= maxIndex; i++) {
        returnValues.push_back(selected[i]);
      }
    } else if (scoreType == "neg") {
      for (int i = k - 1; i >= minIndex; i--) {
        returnValues.push_back(selected[i]);
      }
    }
  }
  
  return(returnValues);
}

// [[Rcpp::export]]
NumericMatrix GSEACpp(List pathway, NumericVector scoresNull, int nGenes, int nReplicate, NumericVector scoreObserved, std::string scoreType, int maxExtreme, int maxPermutation, double p) {
  // gene set enrichment analysis using the paradigm of phentoype permutation,
  // which is better than gene permutation, such as used in fgsea
  //   
  // input
  // pathway: a list of vectors defining gene sets with integers specifying the
  //          index of each gene in the scoreNull matrix or scoreObserved
  // scoresNull: a vector from the matrix of gene*replicate scores based on the 
  //           null hypothesis
  // scoreObserved: observed scores, the gene has the same order as in 
  //           scoresNull
  // scoreType: "pos", "neg", or "std". "pos", "neg" are for one-sided tests. 
  //            "pos" tests whether the pathways have higher enrichment scores, 
  //            "neg" tests whether the pathways have lower enrichment scores,
  //            "std" is the standard used in the original GSEA, tests whether
  //            the pathways have more extreme enrichemnt scores, can be in 
  //            either direction
  
  int nPathway = pathway.length();
  NumericVector pathwayPvalue(nPathway);
  NumericVector nExtreme(nPathway);
  NumericVector nTotal(nPathway);
  LogicalVector finished(nPathway);
  finished.fill(false);
  NumericVector ES(nPathway);
  int nSample = maxPermutation;
  NumericVector scoreSampled(nGenes);
  int nDone = 0;
  
  // calculate the observed ES
  IntegerVector orderIndex = sort_indexes_decreasing(scoreObserved);
  NumericVector scoreObservedSorted = scoreObserved[orderIndex];
  // <debug>
  // cout << nPathway << endl;
  // cout << scoreObservedSorted[0] << " " << scoreObservedSorted[1] << endl;
  // </debug>
  LogicalVector inSet(nGenes);
  LogicalVector inSetReordered(nGenes);
  for (int i = 0; i < nPathway; i++) {
    inSet.fill(false);
    IntegerVector index = pathway[i];
    inSet[index - 1] = true; 
    inSetReordered = inSet[orderIndex];
    ES[i] = calculateESCppFast(scoreObservedSorted, inSetReordered, scoreType = scoreType, p, false)[0];
    
    // <debug>
    // cout << "pathway " << i << endl;
    // for (int j = 0; j < index.size(); j ++ )
    //   cout << " " << index[j];
    // cout << endl;
    // cout << ES[i] << endl;
    // </debug>
  }
  
  // <debug>
  // return 0;
  // </debug>
  
  // calculate the extreme ES in the null
  for (int j = 0; j < nSample; j++) {
    if ((j + 1) % 1000 == 0) {
      cout << j + 1 << endl;
    }
    // sample one sample of scores from all simulated null
    IntegerVector indexOfReplicate = seq(1, nReplicate);
    
    IntegerVector sampledReplicate = RcppArmadillo::sample(indexOfReplicate, nGenes, true);
    // <debug>
    // IntegerVector sampledReplicate = rep(1, nGenes);
    // </debug>
    
    IntegerVector index = (sampledReplicate - 1) * nGenes + seq(1, nGenes) - 1;
    NumericVector scoreSampled = scoresNull[index];
    IntegerVector orderIndex = sort_indexes_decreasing(scoreSampled);
    // <debug>
    // IntegerVector orderIndex = seq(1, nGenes) - 1;
    // </debug>
    NumericVector scoreSampledSorted = scoreSampled[orderIndex];
    
    // <debug>
    // NumericVector scoreSampledSorted = scoresNull[seq(1, nGenes)];
    // </debug>
    
    for (int i = 0; i < nPathway; i++) {
      // if (i % 100 == 0) {
      //   cout << i << endl;
      // }
      if (!finished[i]) {
        inSet.fill(false);
        IntegerVector index = pathway[i];
        inSet[index - 1] = true;
        inSetReordered = inSet[orderIndex];
        double ESNull = calculateESCppFast(scoreSampledSorted, inSetReordered, scoreType = scoreType, p, false)[0];
        
        // <debug>
        // double ESNull = 0;
        // </debug>
        
        if (scoreType == "pos") {
          if (ESNull >= ES[i]) {
            nExtreme[i] = nExtreme[i] + 1;
          }
          nTotal[i] = nTotal[i] + 1;
        } else if (scoreType == "neg") {
          if (ESNull <= ES[i]) {
            nExtreme[i] = nExtreme[i] + 1;
          }
          nTotal[i] = nTotal[i] + 1;
        } else if (scoreType == "std") {
          if (ES[i] >= 0) {
            if (ESNull >= ES[i] && ESNull >= 0) {
              nExtreme[i] = nExtreme[i] + 1;
            } 
            if (ESNull >= 0) {
              nTotal[i] = nTotal[i] + 1;
            }
          } else if (ES[i] < 0) {
            if (ESNull <= ES[i] && ESNull < 0) {
              nExtreme[i] = nExtreme[i] + 1;
            }
            if (ESNull < 0) {
              nTotal[i] = nTotal[i] + 1;
            }
          }
        }
        if ( nExtreme[i] >= maxExtreme ) {
          finished[i] = true;
          nDone = nDone + 1;
        }
      }
    }
    
    if (nDone == nPathway) {
      break;
    }
  }
  
  for (int i = 0; i < nPathway; i++) {
    pathwayPvalue[i] = nExtreme[i] / nTotal[i];
  }
  
  // // <debug>
  // cout << "before the last statement" << endl;
  // // </debug>
  NumericMatrix  pathwayResult = cbind(ES, pathwayPvalue, nExtreme, nTotal);
  
  return pathwayResult;
}


