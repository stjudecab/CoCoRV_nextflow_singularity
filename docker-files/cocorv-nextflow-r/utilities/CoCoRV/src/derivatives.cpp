#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector calculateProbCpp(NumericMatrix p) {
// calculate the probability of no alternate alleles, 1 alternate alleles, 
// 2 alternate alleles and two heterogyzous genotypes
  
// input
// p: either a 1*nSNP matrix of AFs or a 3 * nSNP genotype frequencies
// frequencies.
// When p is a 3*nSNP matrix
//   each row corresponds to genotype RR, RA, AA
//   each column is a variant
  
  if (!((p.nrow() == 1 || p.nrow() == 3) && p.ncol() >= 1)) {
    std::cout << "wrong dimension of the input frequencies" << std::endl;
    exit(1);
  }
  
  int nSNP = p.ncol();
  NumericMatrix geno(3, nSNP);
  if (p.nrow() == 1) {
    for (int j = 0; j < nSNP; j++) {
      geno(0, j) = std::pow(1 - p(0, j), 2);
      geno(1, j) = 2 * (1 - p(0, j)) * p(0, j);
      geno(2, j) = p(0, j) * p(0, j);
    }
  } else {
    geno = p;
  }
  
  double p0, p1, p2, pTwoHets, p2h, pCHets;
  if (nSNP == 1) {
    p0 = geno(0, 0); 
    p1 = 1 - geno(0, 0); 
    p2 = geno(2, 0);
    pTwoHets = 0;
    p2h = p2;
    pCHets = 0;
  } else {
    p0 = 0;
    for (int j = 0; j < nSNP; j++) {
      p0 += std::log(geno(0, j));
    }
    p0 = std::exp(p0);
    
    p1 = 0;
    for (int j = 0; j < nSNP; j++) {
      double prod = std::log(1 - geno(0, j));
      for (int i = 0; i < nSNP; i++) {
        if (i != j) {
          prod += std::log(geno(0, i));
        }
      }
      prod = std::exp(prod);
      p1 += prod;
    }
    
    // exactly one homozyous alternate
    double pOneHomAlt = 0;
    for (int j = 0; j < nSNP; j++) {
      double prod = std::log(geno(2, j));
      for (int i = 0; i < nSNP; i++) {
        if (i != j) {
          prod += std::log(geno(0, i));
        }
      }
      prod = std::exp(prod);
      pOneHomAlt += prod;
    }
    
    // two heterozygous 
    pTwoHets = 0;
    if (nSNP == 2) {
      pTwoHets = (1 - geno(0, 0)) * (1 - geno(0, 1));
    } else {
      for (int j = 0; j < nSNP - 1; j++) {
        for (int i = j + 1; i < nSNP; i++) {
          double prod = std::log(1 - geno(0, j)) + std::log(1 - geno(0, i));
          for (int k = 0; k < nSNP; k++) {
            if (k != i && k != j) {
              prod += std::log(geno(0, k));
            }
          }
          prod = std::exp(prod);
          pTwoHets += prod;
        }
      }
    }
    
    p2 = pOneHomAlt + pTwoHets;
    pCHets = pTwoHets / 2;
    p2h = pOneHomAlt + pCHets;
  }
  
  NumericVector out = NumericVector::create(p0, p1, p2, pTwoHets, p2h, pCHets);
  return(out);
}

// [[Rcpp::export]]
NumericVector oddsRatioToParCpp(NumericVector ORPar) {
// convert the odds ratio based parameters to the probability parameters
  
// input
// ORPar: log10(AF), the log10 ratio of AF between the second and and 
//  the first, log10 odds ratio
  
  double s = pow(10, ORPar[0]);
  double t = pow(10, ORPar[1]);
  double theta = pow(10, ORPar[2]);
  
  double h = (s + t) * ((s + t) * (1 - theta) - 2) + 4 * s * t * theta;
  // double l = pow((s + t) * (1 - theta) - 1, 2) + 
  //       4 * (1 - theta) * s * t * theta;
  double l = 1.0 + (1 - theta) * (pow(s + t, 2) - theta * pow(s - t, 2) - \
        2 * (s + t));
  double p11 = 0.5 * (s + t) + 0.5 * h / (sqrt(l) + 1);
  p11 = std::max(0., p11);
  
  double p10 = std::max(0., s - p11);
  double p01 = std::max(0., t - p11);
  
  NumericVector out = NumericVector::create(p11, p10, p01);
  
  return(out);
}

// [[Rcpp::export]]
double dmultinomCpp(IntegerVector observed, NumericVector probability,
                           bool useLog = false) {
// make sure all probabilities sum to 1 

  int n = std::accumulate(observed.begin(), observed.end(), 0);
  double v = std::lgamma(n + 1);
  for (int i = 0; i < observed.size(); i++) {
    v = v - std::lgamma(observed[i] + 1) + \
         observed[i] * std::log(probability[i]);
  }
  
  if (!useLog) {
    v = std::exp(v);
  }
  
  return(v);
}

// [[Rcpp::export]]
double negativeLogLikCpp(NumericVector par, IntegerMatrix data) {
// the log likelihood of the observed data x, y, n

// input
// ORPar: the three parameters including the odds ratio, see oddsRatioToPar

  IntegerVector x = data(_, 0);
  IntegerVector y = data(_, 1);
  IntegerVector n = data(_, 2);
  int nGroup = x.size();
  double logL = 0;
  NumericVector prob = oddsRatioToParCpp(par);
  double p00 = 1 - prob[0] - prob[1] - prob[2];
  p00 = std::max(0., p00);
  NumericVector probability = NumericVector::create(prob[0], \
                          prob[1], prob[2], p00);
  IntegerVector observed(4);

  for (int i = 0; i < nGroup; i++) {
    int rMin = std::max(0, x[i] + y[i] - n[i]);
    int rMax = std::min(x[i], y[i]);
    if (rMin > rMax) {
      std::cout << "not compatible summary counts" << std::endl;
      return(NumericVector::get_na());
    }
    int rCount = rMax - rMin + 1;
    NumericVector sumPVector(rCount);
    int j = 0;
    for (int r = rMin; r <= rMax; r++) {
      double eps = 0;
      observed[0] = r + eps;
      observed[1] = x[i] - r + eps;
      observed[2] = y[i] - r + eps;
      observed[3] = n[i] - y[i] - x[i] + r + eps;
      sumPVector[j] = dmultinomCpp(observed, probability, true);
      j++;
    }
    // robust calculation of log(sumP): sumP is in still log scale
    double logPMax = max(sumPVector);
    double logSumP = log(sum(exp(sumPVector - logPMax))) + logPMax;
    logL = logL + logSumP;
  }

  return(-logL);
}

// [[Rcpp::export]]
double negativeLogLikGenotypeCpp(NumericVector par, NumericVector data) {
  // calculate the -log likelihood of diploid data
  
  double logL = 0;
  NumericVector prob = oddsRatioToParCpp(par);
  
  double pab = prob[0];
  double paB = prob[1];
  double pAb = prob[2];
  double pAB = std::max(0., 1 - prob[0] - prob[1] - prob[2]);
  
  double logpab = 0;
  double logpaB = 0;
  double logpAb = 0;
  double logpAB = 0;
  if (pab > 0) {
    logpab = std::log(pab);
  }
  if (paB > 0) {
    logpaB = std::log(paB);
  }
  if (pAb > 0) {
    logpAb = std::log(pAb);
  }
  if (pAB > 0) {
    logpAB = std::log(pAB);
  }
  
  NumericVector logf = NumericVector::create(logpAB * 2, \
                                             std::log(2) + logpAB + logpAb, \
                                             logpAb * 2, \
                                             std::log(2) + logpAB + logpaB, \
                      std::log(2) + std::log(pAB * pab + pAb * paB), \
                      std::log(2) + logpAb + logpab, \
                      logpaB * 2, std::log(2) + logpaB + logpab, \
                      logpab * 2);
  
  for (int i = 0; i < 9; i++) {
    if (data[i] > 0) {
      logL = logL + data[i] * logf[i];
    }
  }
  
  return(-logL);
}

NumericVector h2stthetaCpp(NumericVector x) {
  double s = x[0];
  double t = x[1];
  double theta = x[2];
  
  double tos = 2 * (s + t) * (1 - theta) + 4 * t * theta - 2;
  double tot = 2 * (s + t) * (1 - theta) + 4 * s * theta - 2;
  double totheta = -(s + t) * (s + t) + 4 * s * t;
  
  return(NumericVector::create(tos, tot, totheta));
}

NumericVector l2stthetaCpp(NumericVector x) {
  double s = x[0];
  double t = x[1];
  double theta = x[2];
  
  double tos = 2 * ((s + t) * (1 - theta) - 1) * (1 - theta) + 
    4 * (1 - theta) * t * theta;
  double tot = 2 * ((s + t) * (1 - theta) - 1) * (1 - theta) + 
    4 * (1 - theta) * s * theta;
  double totheta = -2 * ((s + t) * (1 - theta) - 1) * (s + t) + 4 * s * t - 
    8 * s * t * theta;
  
  return(NumericVector::create(tos, tot, totheta));
}

NumericVector a2stthetaCpp(NumericVector x) {
  double s = x[0];
  double t = x[1];
  double theta = x[2];
  
  double h = (s + t) * ((s + t) * (1 - theta) - 2) + 4 * s * t * theta;
  double l = std::pow((s + t) * (1 - theta) - 1, 2) + \
             4 * (1 - theta) * s * t * theta;
  NumericVector hgrad = h2stthetaCpp(x);
  NumericVector lgrad = l2stthetaCpp(x);
  NumericVector agrad(3);
  NumericVector part1 = NumericVector::create(0.5, 0.5, 0);
  for (int i = 0; i < 3; i++) {
    agrad[i] = part1[i] + 0.5 * (hgrad[i] * (std::sqrt(l) + 1) - 
    0.5 * h / sqrt(l) * lgrad[i]) / (l + 2 * std::sqrt(l) + 1);
  }
  return(agrad);
}

// [[Rcpp::export]]
NumericVector negativeLogLikGradCpp(NumericVector par, IntegerMatrix data) {
// the gradient of the log likelihood of the observed data x, y, n
  
// input
// ORPar: the three parameters including the odds ratio, see oddsRatioToPar
  
  IntegerVector x = data(_, 0);
  IntegerVector y = data(_, 1);
  IntegerVector n = data(_, 2);
  int nGroup = x.size();
  NumericVector logLGrad(3);
  logLGrad[0] = 0;
  logLGrad[1] = 0;
  logLGrad[2] = 0;
  NumericVector prob = oddsRatioToParCpp(par);
    
  double s = pow(10, par[0]);
  double t = pow(10, par[1]);
  double theta = pow(10, par[2]);
  NumericVector par2 = NumericVector::create(s, t, theta);
    
  double a = prob[0];
  double b = prob[1];
  double c = prob[2];
  double d = std::max(0., 1 - a - b - c);
    
  NumericVector probability = NumericVector::create(a, b, c, d);
    
  NumericVector aDerivative = a2stthetaCpp(par2);
  NumericVector bDerivative(3);
  bDerivative[0] =  1 - aDerivative[0];
  bDerivative[1] =  - aDerivative[1];
  bDerivative[2] =  - aDerivative[2];
  NumericVector cDerivative(3);
  cDerivative[0] =  - aDerivative[0];
  cDerivative[1] =  1 - aDerivative[1];
  cDerivative[2] =  - aDerivative[2];
  NumericVector dDerivative(3);
  dDerivative[0] =  -1 + aDerivative[0];
  dDerivative[1] =  -1 + aDerivative[1];
  dDerivative[2] =  aDerivative[2];
    
  double loga = 0;
  double logb = 0;
  double logc = 0;
  double logd = 0;
  if (a > 0) {
    loga = std::log(a);
  }
  if (b > 0) {
    logb = std::log(b);
  }
  if (c > 0) {
    logc = std::log(c);
  }
  if (d > 0) {
    logd = std::log(d);
  }
  
  IntegerVector observed(4);
  for (int i = 0; i < nGroup; i++) {
    int rMin = std::max(0, x[i] + y[i] - n[i]);
    int rMax = std::min(x[i], y[i]);
    if (rMin > rMax) {
      std::cout << "not compatible summary counts" << std::endl;
      return(NumericVector::get_na());
    }
    double sumP = 0;
    NumericVector gradPerRow(3);
    gradPerRow[0] = 0;
    gradPerRow[1] = 0;
    gradPerRow[2] = 0;
    double logkInit = std::lgamma(n[i] + 1);
    for (int r = rMin; r <= rMax; r++) {
      // calculate the sum of probability
      double eps = 0;
      observed[0] = r + eps;
      observed[1] = x[i] - r + eps;
      observed[2] = y[i] - r + eps;
      observed[3] = n[i] - y[i] - x[i] + r + eps;
      double fullLoglik = dmultinomCpp(observed, prob = probability, false);
      sumP = sumP + fullLoglik;
      
      // calculate the derivatives
      int o1 = observed[0];
      int o2 = observed[1];
      int o3 = observed[2];
      int o4 = observed[3];
      double logk = logkInit;
      for (int j = 0; j < 4; j++) {
        logk = logk - std::lgamma(observed[j] + 1);
      }
      double part1 = 0;
      double part2 = 0;
      double part3 = 0;
      double part4 = 0;

      for (int j = 0; j < 2; j++) {      
        if (o1 > 0) {
          part1 = std::log(o1) + loga * (o1 - 1) + o2 * logb + \
               o3 * logc + o4 * logd + logk;
          part1 = exp(part1);
        }
        if (o2 > 0) {
          part2 = loga * o1 + (o2 - 1) * logb + std::log(o2) + o3 * logc + \
            o4 * logd + logk;
          part2 = exp(part2);
        }
        if (o3 > 0) {
          part3 = loga * o1 + o2 * logb + std::log(o3) + (o3 - 1) * logc + \
              o4 * logd + logk;
          part3 = exp(part3);
        }
        if (o4 > 0) {
          part4 = loga * o1 + o2 * logb + o3 * logc + (o4 - 1) * logd + \
              std::log(o4) + logk;
          part4 = exp(part4);
        }
        for (j = 0; j < 3; j++) {
          gradPerRow[j] = gradPerRow[j] + (part1 * aDerivative[j] + \
                    part2 * bDerivative[j] + \
                    part3 * cDerivative[j] + part4 * dDerivative[j]) * \
                    par2[j] * std::log(10);
        }
      }
    }
    if (sumP > 0) {
      for (int j = 0; j < 3; j++) {
        gradPerRow[j] = gradPerRow[j] / sumP;
        logLGrad[j] = logLGrad[j] + gradPerRow[j];
      }
    }
  }
  
  return(logLGrad * (-1));
}

// [[Rcpp::export]]
NumericVector negativeLogLikGenotypeGradCpp(NumericVector par, NumericVector data) {
  // the gradient of the log likelihood of the observed data x, y, n

  // input
  // ORPar: the three parameters including the odds ratio, see oddsRatioToPar

  NumericVector logLGrad(3);
  logLGrad[0] = 0;
  logLGrad[1] = 0;
  logLGrad[2] = 0;
  NumericVector prob = oddsRatioToParCpp(par);

  double s = pow(10, par[0]);
  double t = pow(10, par[1]);
  double theta = pow(10, par[2]);
  NumericVector par2 = NumericVector::create(s, t, theta);

  double a = prob[0];
  double b = prob[1];
  double c = prob[2];
  double d = std::max(0., 1 - a - b - c);

  NumericVector probability = NumericVector::create(a, b, c, d);

  NumericVector aDerivative = a2stthetaCpp(par2);
  NumericVector bDerivative(3);
  bDerivative[0] =  1 - aDerivative[0];
  bDerivative[1] =  - aDerivative[1];
  bDerivative[2] =  - aDerivative[2];
  NumericVector cDerivative(3);
  cDerivative[0] =  - aDerivative[0];
  cDerivative[1] =  1 - aDerivative[1];
  cDerivative[2] =  - aDerivative[2];
  NumericVector dDerivative(3);
  dDerivative[0] =  -1 + aDerivative[0];
  dDerivative[1] =  -1 + aDerivative[1];
  dDerivative[2] =  aDerivative[2];

  NumericVector f = NumericVector::create(d * d, 2 * c * d, c * c, \
                          2 * b * d, 2 * (a * d + b * c), 2 * a * c, \
                          b * b, 2 * a * b, a * a);
  NumericVector f2a = NumericVector::create(0, 0, 0, 0, 2 * d, 2 * c, \
                                            0, 2 * b, 2 * a);
  NumericVector f2b = NumericVector::create(0, 0, 0, 2 * d, 2 * c, 0, \
                                              2 * b, 2 * a, 0);
  NumericVector f2c = NumericVector::create(0, 2 * d, 2 * c, 0, 2 * b, 2 * a, \
                                            0, 0, 0);
  NumericVector f2d = NumericVector::create(2 * d, 2 * c, 0, 2 * b, 2 * a, 0, \
                                              0, 0, 0);

  for (int j = 0; j < 3; j++) {
    for (int i = 0; i < 9; i++) {
      if (data[i] > 0) {
        logLGrad[j] = logLGrad[j] +  \
          data[i] / f[i] * (f2a[i] * aDerivative[j] + f2b[i] * bDerivative[j] \
                      + f2c[i] * cDerivative[j] + f2d[i] * dDerivative[j]);
      }
    }
    logLGrad[j] = logLGrad[j] * par2[j] * std::log(10);
  }

  return(logLGrad * (-1));
}