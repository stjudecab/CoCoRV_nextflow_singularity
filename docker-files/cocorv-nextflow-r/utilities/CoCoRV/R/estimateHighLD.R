negativeLogLikGenotypeNull = function(par, data, log10OR = 0) {
  value = negativeLogLikGenotypeCpp(c(par, log10OR), data)
  return(value)
}

negativeLogLikGenotypeGradNull = function(par, data, log10OR = 0) {
  value = negativeLogLikGenotypeGradCpp(c(par, log10OR), data)
  return(value[1:2])
}

negativeLogLikNull = function(par, data, log10OR = 0) {
  value = negativeLogLikCpp(c(par, log10OR), data)
  return(value)
}

negativeLogLikGradNull = function(par, data, log10OR = 0) {
  value = negativeLogLikGradCpp(c(par, log10OR), data)
  return(value[1:2])
}

fitNull = function(data, log10OR = 0, genotype = F) {
  x = data[, 1]
  y = data[, 2]
  if (genotype) {
    n = rep(2, length(x))
  } else {
    n = data[, 3]
  }
  AF1 = min(max(c(1e-5, sum(x) / sum(n))), 1 - 1e-5)
  AF2 = min(max(c(1e-5, sum(y) / sum(n))), 1 - 1e-5)
  
  result = NULL
  parInit = c(log10(AF1), log10(AF2))
  if (genotype) {
    opts = list(
      algorithm = "NLOPT_LD_MMA",
      xtol_rel = 1e-04)
    n1 = sum(x == 0 & y == 0)
    n2 = sum(x == 0 & y == 1)
    n3 = sum(x == 0 & y == 2)
    n4 = sum(x == 1 & y == 0)
    n5 = sum(x == 1 & y == 1)
    n6 = sum(x == 1 & y == 2)
    n7 = sum(x == 2 & y == 0)
    n8 = sum(x == 2 & y == 1)
    n9 = sum(x == 2 & y == 2)
    genotypeData = c(n1, n2, n3, n4, n5, n6, n7, n8, n9)
    result = nloptr(x0 = parInit, 
                         eval_f = negativeLogLikGenotypeNull,
                         eval_grad_f = negativeLogLikGenotypeGradNull,
                         lb = c(-Inf, -Inf), 
                         ub = c(0, 0),
                         opts = opts,
                         log10OR = log10OR,
                         data = genotypeData)
  } else {
    opts = list(
      algorithm = "NLOPT_LD_SLSQP", 
      xtol_rel = 1e-04)
    result = nloptr(x0 = parInit, 
                         eval_f = negativeLogLikNull,
                         eval_grad_f = negativeLogLikGradNull,
                         lb = c(-Inf, -Inf), 
                         ub = c(0, 0),
                         opts = opts,
                         log10OR = log10OR,
                         data = data)
  }
  
  solution = result$solution
  logLik = -result$objective
  
  return(list(x = solution, logLik = logLik))
}

oddsRatioMLE = function(data, ORUpper = 10^10, ORLower = 0, genotype = F) {
  # MLE of OR
  # assume ORUpper is at least 1
  
  # input
  # data:
  # ORUpper: the maximal OR threshold
  # ORLower: the minimal OR threshold
  
  stopifnot(ORUpper >= 1)

  x = data[, 1]
  y = data[, 2]
  if (genotype) {
    n = rep(2, length(x))
  } else {
    n = data[, 3]
  }
  AF1 = min(max(c(1e-5, sum(x) / sum(n))), 1 - 1e-5)
  AF2 = min(max(c(1e-5, sum(y) / sum(n))), 1 - 1e-5)
  
  # convert to log10 scale
  ORLower = log10(ORLower)
  ORUpper = log10(ORUpper)
  
  #browser()
  result = list()
  bestResult = NA
  minNegLikelihood = Inf
  if (ORUpper > 0) {
    if (ORLower <= 3) {
      parInitSet = cbind(log10(AF1), log10(AF2), 
                       seq(from = max(-2, ORLower), by = 1, 
                           to = min(ORUpper, 3)))
    } else {
      parInitSet = cbind(log10(AF1), log10(AF2), 
                         c(ORLower, ORLower + 1))
    }
  } else {
    parInitSet = cbind(log10(AF1), log10(AF2), 
                       seq(from = max(-5, ORLower), by = 1, 
                           to = 0))
  }
  for (i in 1:dim(parInitSet)[1]) {
    parInit = parInitSet[i, ]
    if (genotype) {
      opts = list(
        #algorithm = "NLOPT_LD_SLSQP", 
        algorithm = "NLOPT_LD_MMA",
        #algorithm = "NLOPT_LN_BOBYQA",
        xtol_rel = 1e-04)
      n1 = sum(x == 0 & y == 0)
      n2 = sum(x == 0 & y == 1)
      n3 = sum(x == 0 & y == 2)
      n4 = sum(x == 1 & y == 0)
      n5 = sum(x == 1 & y == 1)
      n6 = sum(x == 1 & y == 2)
      n7 = sum(x == 2 & y == 0)
      n8 = sum(x == 2 & y == 1)
      n9 = sum(x == 2 & y == 2)
      genotypeData = c(n1, n2, n3, n4, n5, n6, n7, n8, n9)
      result[[i]] = nloptr(x0 = parInit, 
                           eval_f = negativeLogLikGenotypeCpp,
                           eval_grad_f = negativeLogLikGenotypeGradCpp,
                           # eval_f = negativeLogLikGenotype,
                           #lb = rep(0, 3), ub = rep(1, 3),
                           lb = c(-Inf, -Inf, ORLower), 
                           ub = c(0, 0, ORUpper),
                           #lb = c(-Inf, -Inf, -Inf), ub = c(0, 0, 0),
                           opts = opts,
                           data = genotypeData)
    } else {
      opts = list(
        algorithm = "NLOPT_LD_SLSQP", 
        #algorithm = "NLOPT_LD_MMA",
        #algorithm = "NLOPT_LN_BOBYQA",
        xtol_rel = 1e-04)
      result[[i]] = nloptr(x0 = parInit, 
                         eval_f = negativeLogLikCpp,
                         eval_grad_f = negativeLogLikGradCpp,
                         #eval_f = negativeLogLik,
                         #eval_grad_f = negativeLogLikGrad,
                         #lb = rep(0, 3), ub = rep(1, 3),
                         lb = c(-Inf, -Inf, ORLower), 
                         ub = c(0, 0, ORUpper),
                         #lb = c(-Inf, -Inf, -Inf), ub = c(0, 0, 0),
                         opts = opts,
                         data = data)
      if (is.na(result[[i]]$objective)) {
        opts = list(
          algorithm = "NLOPT_LD_MMA",
          xtol_rel = 1e-04) 
        result[[i]] = nloptr(x0 = parInit, 
                             eval_f = negativeLogLikCpp,
                             eval_grad_f = negativeLogLikGradCpp,
                             #eval_f = negativeLogLik,
                             #eval_grad_f = negativeLogLikGrad,
                             #lb = rep(0, 3), ub = rep(1, 3),
                             lb = c(-Inf, -Inf, ORLower), 
                             ub = c(0, 0, ORUpper),
                             #lb = c(-Inf, -Inf, -Inf), ub = c(0, 0, 0),
                             opts = opts,
                             data = data)
      }
      if (is.na(result[[i]]$objective)) {
        opts = list(
          algorithm = "NLOPT_LN_BOBYQA",
          xtol_rel = 1e-04) 
        result[[i]] = nloptr(x0 = parInit, 
                             eval_f = negativeLogLikCpp,
                             eval_grad_f = negativeLogLikGradCpp,
                             #eval_f = negativeLogLik,
                             #eval_grad_f = negativeLogLikGrad,
                             #lb = rep(0, 3), ub = rep(1, 3),
                             lb = c(-Inf, -Inf, ORLower), 
                             ub = c(0, 0, ORUpper),
                             #lb = c(-Inf, -Inf, -Inf), ub = c(0, 0, 0),
                             opts = opts,
                             data = data)
      }
    }
    if (is.na(bestResult) || !is.na(result[[i]]$objective) & 
        result[[i]]$objective < minNegLikelihood) {
      bestResult = i
      minNegLikelihood = result[[i]]$objective
      status = result[[i]]$status
    }
  }
  solution = result[[bestResult]]$solution
  
  return(list(x = solution, logLik = -minNegLikelihood, status = status))
}

#' This function uses a set of summary counts to estimate/test the linkage 
#' disequilibrium
#' 
#' It can be applied to both the summary counts or full genotypes, where the
#' latter is a special case of pooling two haplotypes together

#' @param data a matrix of 3 or 2 columns. The first two columns are the pooled
#' alternate allele counts of the two variants. The third column is the number
#' of haplotypes in total. If genotype data is used, the third column can be 
#' omitted
#' @param ORThreshold the threshold used in the test
#' @param alternative either "greater" or "less". If "greater", the alternative
#' hypothesis is that the odds ratio theta > ORThreshold. If "less", then it is
#' theta < ORThreshold
#' @param genotype If true, then assume the data is generated using full 
#' genotype data, i.e., n = 2 for all samples. Note that in this case, each row
#' corresponds to an individual
#' @param ACThreshold the minimal alternate allele count for odds ratio 
#' estimation or test

#' @return A vector of 12 values: pvalue of the 
#' likelihood ratio test and the estimated odds ratio, r2, D', log likelihood 
#' of the full and the null model, the next four are estimated 
#' probalities/frequencies of the 4 haplotypes: p11, p10, p01, p00, and the 
#' return status of the optimization for the full and the null model, where
#' 1-4 indicate success, other values indicate failing to converge  

#' @export
LDTest = function(data, ORThreshold = 1,
                  alternative = "greater",
                  genotype = F,
                  ACThreshold = 1) {
  if (genotype) {
    # check the number of total ACs to avoid too few counts
    totalAC = sum(data[, 1:2], na.rm = T)
    if (totalAC < ACThreshold) {
      return(c(1, rep(NA, 7)))
    }
  }
  
  # impute 0 to all NAs in the data, this is likely for genotype data
  data[is.na(data)] = 0
  
  if (sum(data[, 1:2]) == 0) {
    print("LD test not performed, need >=1 alternate alleles")
    return(rep(NA, 8))
  }
  
  #browser()
  # full model
  result1 = oddsRatioMLE(data, ORUpper = Inf, genotype = genotype)
  solution = result1$x
  
  # direct one-sided test 
  result0 = NA
  if (alternative == "greater") {
    # restrict to the null theta <= ORThreshold
    result0 = oddsRatioMLE(data, ORUpper = ORThreshold, genotype = genotype)
  } else if (alternative == "less") {
    result0 = oddsRatioMLE(data, ORLower = ORThreshold, genotype = genotype)
  } else if (alternative == "two.sided" ) {
    result0 = oddsRatioMLE(data, ORUpper = ORThreshold, 
                           ORLower = ORThreshold, genotype = genotype)
  } else {
    stop("wrong input parameter for alternative, must be greater, less, or
         two.sided")
  }
  #
  # LRT
  zLRT = 2 * (result1$logLik - result0$logLik)
  pLRT = pchisq(zLRT, df = 1, lower.tail = F)
  
  # signed root of LR
  # calculate the logLik of the null
  # L0 = fitNull(data, log10OR = log10(ORThreshold), genotype = genotype)
  # LRT = 2 * (result1$logLik - L0$logLik)
  # if (LRT < 0) {
  #   LRT = 0
  # }
  # w = sign(solution[3] - log10(ORThreshold)) * sqrt(LRT)
  # pLRT = NA
  # if (alternative == "greater") {
  #   pLRT = pnorm(w, 0, 1, lower.tail = F)
  # } else if (alternative == "less") {
  #   pLRT = pnorm(w, 0, 1, lower.tail = T)
  # } else {
  #   stop("wrong input parameter for alternative")
  # }

  #browser()
  # calculate the 95% confidence interval
  # lowerCI = NA
  # upperCI = NA
  # library(numDeriv)
  # hessianMatrix = hessian(
  #                         #func = negativeLogLikBasedOnEM, 
  #                         func = negativeLogLikCpp, 
  #                         x = solution, data = data)
  # eigenResult = eigen(hessianMatrix)
  # lambda = eigenResult$values
  # V = eigenResult$vectors
  # lambdaInverse = 1 / lambda
  # # information matrix of logL as the hessian of -logL
  # if (min(lambdaInverse) < 0) {
  #   if (min(lambda) < -0.1) {
  #     print("The smallest eigen value of the information matrix < -0.1")
  #   }
  #   #browser()
  #   lambdaInverse[lambdaInverse < 0] = -lambdaInverse[lambdaInverse < 0]
  # }
  # covarianceMatrix = V %*% diag(lambdaInverse) %*% t(V)
  # lowerCI = solution[3] - 1.96 * sqrt(covarianceMatrix[3, 3]) 
  # upperCI = solution[3] + 1.96 * sqrt(covarianceMatrix[3, 3])
  
  log10OR = solution[3]
  OR = 10^log10OR
  
  prob = oddsRatioToParCpp(solution)
  
  r2AndDPrime = r2DPrime(c(prob, 1 - sum(prob))) 
  optimStatusFull = result1$status
  optimStatusNull = result0$status
  result = c(pLRT, OR, r2AndDPrime, result1$logLik, result0$logLik, 
    prob, 1 - sum(prob), optimStatusFull, optimStatusNull)
  names(result) = c("pvalue", "OR", "r2", "DPrime", "logLikFull", "logLikNull", "p11", "p10", "p01", "p00", "optimStatusFull", "optimStatusNull")
  return(result)
}

r2DPrime = function(probabilities) {
  # calculate the r2 and DPrime from the estimated haplotype probabilites
  # assume biallelic for both variants

  # input
  # probabilities: a vector p11, p10, p01 and p00

  # output
  # a vector of r2 and D'

  p11 = probabilities[1]
  p10 = probabilities[2]
  p01 = probabilities[3]
  p00 = probabilities[4]

  p1 = p11 + p10
  q1 = p11 + p01
  D11 = p11 - p1 * q1
  D11Prime = NA
  if (D11 < 0) {
    D11Prime = D11 / max(- p1 * q1, - (1 - p1) * (1 - q1))
  } else {
    D11Prime = D11 / min(p1 * (1 - q1), (1 - p1) * q1)
  }
  
  r2 = D11^2 / (p1 * (1 - p1) * q1 * (1 - q1))
  result = c(r2, D11Prime)
  names(result) = c("r2", "DPrime")
  return(result)  
}