#' This function calculates the re-sampling based FDR 

#' It is from Yekutieli & Benjamini, 1999, but with a minor modification 
#' on the estimated v(p) under the null, instead of m * p, it is directly 
#' estimated from the simulated null pvalues. The upper limit RBH only 
#' considers x = p in Equation 10 of Yekutieli & Benjamini, 1999. 
#' These modifications are also used in the R package FDR-AME 
#' implementing these methods, however, they use < instead of <= when defining
#' r, different from the paper. Here we use <= to define r as in the paper. 

#' @param observed the observed p value
#' @param simulated a matrix of mTest * nPermutation of p values under the null, 
#'   make sure the rows are consistent between the observed p values and the
#'   simulated p-values
#' @param alpha set the 1 - alpha quantile used in the algorithm 

#' @return a matrix of two columns, the first is the point estimate FDR, 
#' the second is the upper limit estimate FDR.

#' @export
RBH = function(observed, simulated, alpha = 0.05) {
  nGene = length(observed)
  stopifnot(nGene > 1)
  
  # calculate r for each p-value
  r = numeric(nGene)
  for (i in 1:nGene) {
    r[i] = sum(observed <= observed[i])
  }
  
  RBH = matrix(NA, nGene, 2)
  colnames(RBH) = c("RBH_point", "RBH_upper_limit")
  
  # sort simulated
  rMatrix = matrix(NA, nGene, dim(simulated)[2])
  for (i in 1:dim(simulated)[2]) {
    simulated[, i] = sort(simulated[, i])
  }
  for (i in 1:dim(simulated)[2]) {
    rMatrix[, i] = findInterval(observed, simulated[, i])
  }
  
  for (i in 1:nGene) {
    p = observed[i]
    # r0 = apply(simulated <= p, MARGIN = 2, FUN = sum)
    r0 = rMatrix[i, ]
    
    r0Quantile = quantile(r0, 1 - alpha)
    r0Mean = mean(r0)
    
    # point estimate
    s = r[i] - r0Mean
    if ( s >= r0Quantile) {
      t = ifelse(r0 + s != 0, r0 / (r0 + s), 0)
      RBH[i, 1] = mean(t)
    } else {
      RBH[i, 1] = mean(r0 > 0)
    }
    
    # upper limit estimate
    s = r[i] - r0Quantile
    if (s > 0) {
      RBH[i, 2] = mean(r0 / (r0 + s))
    } else {
      RBH[i, 2] = mean(r0 > 0)
    }
  }

  # use cummin to enforce non-decreasing fdr for sorted observed p
  o = order(observed, decreasing = T)
  ro = order(o)

  raw = RBH[, 1]
  revised = pmin(1, cummin(raw[o]))[ro]
  RBH[, 1] = revised

  raw = RBH[, 2]
  revised = pmin(1, cummin(raw[o]))[ro]
  RBH[, 2] = revised   
  
  return(RBH)
}

#' calculate RBH based FDR for FET using a 2x2 table or CMH using 2x2xm table
#' If the column number of the input is 4, other count based FDRs for FET is 
#' also calculated

#' @param data a matrix of 4*n 
#' @param alternative should be two.sided, greater, or less as in fisher.test
#' @param nReplication the number of replications in the resampling based 
#' method
#' @param ncore the number of cores to use
#' @param sampling cdf: sampling p-values directly from the cdf of the p-values; 
#'                 table: sampling contingency tables first and then calcualte
#'                 p-values.  

#' @return a matrix of adjusted p-values. If the input is for CMH, then two 
#' methods are returned: RBH_point_estimate, RBH_upper_limit. If the 
#' input is for FET, additoinal FET based FDR are also returned. 

#' @export
countBasedFDR = function(data, alternative = "two.sided",  
                  nReplication = 1e3,
                  ncore = 16, sampling = "cdf") {
  stopifnot(dim(data)[2] %% 4 == 0) # should be 4 * n columns
  if ("data.frame" %in% class(data)) {
    data = as.matrix(data)
  }
  
  # select rows with at least one count for either case or control
  index1 = rowSums(data[, seq(from = 1, by = 2, to = dim(data)[2] - 1)]) > 0
  index2 = rowSums(data[, seq(from = 2, by = 2, to = dim(data)[2])]) > 0
  index = index1 & index2
  if (!all(index)) {
    data = data[index, , drop = F]
  }
  
  if (sampling == "cdf") {
    if (dim(data)[2] > 4) {  ## CMH 
      simuatedResult = simulateNullPvaluesCMH(data, nReplication = nReplication,
                                              alternative = alternative)
    } else {  ## FET 
      # generate the cdf list and sample the p-values directly
      simuatedResult = simulateNullPvaluesFET(data, nReplication = nReplication,
                                              alternative = alternative)
    }
  } else if (sampling == "table") {
    # generate simualted p-values using CoCoRV
    simuatedResult = simulateNullPvalues(data, nReplication = nReplication,
                                         alternative = alternative, 
                                         ncore = ncore)
  }
  observed = simuatedResult$pvalues
  simulated = simuatedResult$nullP
  
  RBHResult = RBH(observed, simulated)
  
  result = cbind(observed, RBHResult)
  colnames(result) = c("pvalue", "RBH_point_estimate", "RBH_upper_limit")
  
  if (dim(data)[2] == 4) { # FET
    # also add other powerful FDR methods
    pvalueAndCDF = fisher.pvalues.support(counts = data, input = "noassoc", 
                                          alternative = alternative)
  } else if (dim(data)[2] %% 4 == 0 && dim(data)[2] > 4) { # CMH
    pvalueAndCDF = cmh.pvalues.support(counts = data, 
                                       alternative = alternative)
  }

  # DBH from package discreteMTP
  DBH = p.discrete.adjust(pvalueAndCDF$raw, pvalueAndCDF$support, 
                                   method = "DBH")
  # step down version of DBH and adaptive DBH from package DiscreteFDR
  DBH.sd = DBH(pvalueAndCDF$raw, pvalueAndCDF$support,
               direction = "sd")$Adjusted
  ADBH.sd = ADBH(pvalueAndCDF$raw, pvalueAndCDF$support,
                 direction = "sd")$Adjusted
  result = cbind(result, DBH, DBH.sd, ADBH.sd)
  
  if (all(index)) {
    return(result)
  } else {
    # fill those filtered with 1
    result2 = matrix(1, length(index), dim(result)[2])
    colnames(result2) = colnames(result)
    result2[index, ] = result
    return(result2)
  }
}


