estimateLambdaEmp = function(pvalues, nullPSortedMean, quantile = 0.95,
                             maxLogP = 30) {
  # estimate lambda based on empirical null P-values  
  
  # input
  # nullPSortedMean: the mean of sorted P-values under the null. 
  
  pvalues[pvalues < 10^(-maxLogP)] = 10^(-maxLogP)
  
  # avoid numeric precision artifact around p value = 1
  pvalues[abs(pvalues - 1) < 1e-12] = 1
  nullPSortedMean[abs(nullPSortedMean - 1) < 1e-12] = 1
  
  observed = sort(-log10(pvalues))
  expected = -log10(nullPSortedMean)
  index95 = 1:round(length(observed) * quantile)
  if (dim(unique(cbind(observed[index95], expected[index95])))[1] <= 20) {
    lambdaEmp = NA
    intercept = NA
  } else {
    fitSummary = summary(lm(observed[index95] ~ expected[index95]))
    if (dim(fitSummary$coefficients)[1] != 2) {
      lambdaEmp = NA
      intercept = NA
    } else {
      lambdaEmp = fitSummary$coefficients[2, 1]
      intercept = fitSummary$coefficients[1, 1]
    }
  }
  
  return(list(lambdaEmp = lambdaEmp, 
              nullPSortedMean = nullPSortedMean,
              intercept = intercept))
}

estimateLambdaQT = function(pvalues, upper = 0.95, lower = 0.5, 
                            removeOne = F, twoPoints = F,
                            maxLogP = 30) {
  ## assume uniform distribution for the expected order statistics 
  ## This implements 3 different ways to estimate lambda
  ## 1) twoPoints: based on the point with the observed pvalue of 1 with the 
  #     highest rank and the pvalue at the 95 percentile among the non-one 
  #     pvalues. Set twoPoints to TRUE to run this method
  ## 2) 0.5-0.95: fit a line to the points within 0.5 and 0.95 quantiles.
  ##     Set removeOne and twoPoints to FALSE to run this method
  ## 3) removeOne: fit a line to the points starting with the point with 
  ##   the observed pvalue of 1 and the highest rank and all points with 
  ##   non-one pvalues until the point with the p-value of 0.95 quantile 
  ##     Set removeOne to TRUE and twoPoints to FALSE to run this method
  
  ## output
  ## b1, a1: b1 is the slope or lambda, a1 is the intercept
  
  pvalues[pvalues < 10^(-maxLogP)] = 10^(-maxLogP)
  pvalues[abs(pvalues - 1) < 1e-15] = 1  # for later == 1 is used
  pvalues = pvalues[order(pvalues, decreasing = T)]
  expected = (length(pvalues):1) / (length(pvalues) + 1)
  
  log10P = -log10(pvalues)
  log10Expected = -log10(expected)
  
  if (twoPoints) {
    y0 = 0
    index1 = which(pvalues == 1)
    stopifnot(length(index1) > 0)
    x0 = log10Expected[max(index1)]
    log10PNonzero = log10P[log10P != 0]
    log10ExpectedNonzero = log10Expected[log10P != 0]
    y1 = log10PNonzero[round(length(log10PNonzero) * upper)]
    x1 = log10ExpectedNonzero[round(length(log10PNonzero) * upper)]
    b1 = (y1 - y0) / (x1 - x0)
    a1 = y0 - b1 * x0
  } else {
    # get the fitted
    if (!removeOne) {
      index = round(length(log10P) * lower):round(
                         length(log10P) * upper)
    } else {
      index1 = which(pvalues == 1)
      if (length(index1) == 0) {
        index = 1:round(length(log10P) * upper)
      } else {
        index = max(index1):round(length(log10P) * upper)
      }
    }
    x = log10Expected[index]
    y = log10P[index]
    coefficients = summary(lm(y ~ x))$coefficients
    a1 = coefficients[1, 1]
    b1 = coefficients[2, 1]
  }
  
  return(c(b1, a1))
}

#' qqplot using the simulated null P values and calculate the FDR designed for
#' discrete counts
#'
#' If the number of columns is 4, it uses the Fisher exact test.
#' If the number of columns if a multiple of 4, it uses the CMH exact test

#' @param data a  matrix or data frame with the number of columns being 4*k:
#'  If k is 1, it uses the Fisher exact test;
#'  if k > 1, it uses the CMH exact test. The columns should be arranged as
#'  follows for each stratified group
#'    counts with rare mutations in case, 
#'    counts without rare mutations in case, 
#'    counts with rare mutations in controls
#'    counts without rare mutations in controls
#' @param quantile the upper quantil of points to use for estimation
#' @param nReplication the number of replications of empirical null P values 
#' @param returnSimulated whether to return simulated null count tables
#' @param returnNullP whether to return simulated null P values
#' @param maxLogP the maximum of -log10(P) to be plotted, larger values will be 
#' set to this maximum value
#' @param nullP the null P values. If not supplied, it will be simulated
#' @param doPlot whether to generate the QQ plot 
#' @param limit set the x and y axis limit  
#' @param seed set the random seed
#' @param plotNULL whether to plot the boxplot of lambda values under the null
#' @param cexRatio set cex values in the plot
#' @param verbose whether verbose or not
#' @param nolambda if true, no lambda will be estimated/plotted 
#' @param alternative see alternative in fisher.test
#' @param ncore the number of cores used
#' @param sampling cdf: sampling p-values directly from the cdf of the p-values; 
#'                 table: sampling contingency tables first and then calculate
#'                 p-values.  
#' @param ... additional parameters to plot


#' @return A list of the following components. 
#' lambdaEmp: the estimated inflation factor of the observed data  
#' lambdaNull: the estimated inflation factors of simulated null data
#' nullPSortedMean: the mean of sorted null P values
#' simulatedCount: the simulated count matrix under the null hypothesis
#' nullP: the matrix of null P-values from the simulated counts
#' pvalues: the observed P-values
#' FDRResult: the FDR results if returnFDR is true
#' indexKept: the index of rows after filtering rows with the total counts 
#'   being 0 
#' intercept: the estimated intercept in the QQ plot

#' @export
qqplotHGAndFDR = function(data, quantile = 0.95, 
                    nReplication = 1000, returnFDR = T,
                    returnSimulated = F, returnNullP = F,
                    maxLogP = 30,
                    nullP = NULL, doPlot = T, limit = NULL, 
                    seed = NULL, plotNULL = F, cexRatio = 1.5, 
                    verbose = T, nolambda = F, alternative = "two.sided", 
                    ncore = 1, sampling = "cdf",
                    ...) {

  # check NAs
  if (any(rowSums(is.na(data)) > 0)) {
    print("NA in the data, please remove them first")
    stop()
  }
  
  # select rows with at least one count for either case or control
  index1 = rowSums(data[, seq(from = 1, by = 2, to = dim(data)[2] - 1)]) > 0
  index2 = rowSums(data[, seq(from = 2, by = 2, to = dim(data)[2])]) > 0
  index = which(index1 & index2)
  if (length(index) < dim(data)[1]) {
    data = data[index, , drop = F]
  }
  
  # simulate data based on the null hypothesis
  simulatedCount = NA
  if (!is.null(nullP)) {
    nReplication = dim(nullP)[2]
    # get the observed P values
    simulatedResult = simulateNullPvalues(data, alternative = alternative,
                                          nReplication = nReplication, 
                                          returnSimulated = F, 
                                          returnNullP = F,
                                          verbose = verbose, seed = seed,
                                          ncore = ncore)
    pvalues = simulatedResult$pvalues  # observed p values
  } else {
    if (!returnSimulated && sampling == "cdf") {
      if (dim(data)[2] == 4) { # FET
        simulatedResult = simulateNullPvaluesFET(data, 
                                            alternative = alternative,
                                            nReplication = nReplication, 
                                            seed = seed)
      } else if (dim(data)[2] %% 4 == 0 && dim(data)[2] > 4) { # CMH
        simulatedResult = simulateNullPvaluesCMH(data, 
                                                 alternative = alternative,
                                                 nReplication = nReplication, 
                                                 seed = seed)
      }
      nullP = simulatedResult$nullP
      pvalues = simulatedResult$pvalues
    } else {
      simulatedResult = simulateNullPvalues(data, alternative = alternative,
                        nReplication = nReplication, 
                        returnSimulated = returnSimulated, 
                        returnNullP = T,
                        verbose = verbose, seed = seed, ncore = ncore)
      nullP = simulatedResult$nullP
      pvalues = simulatedResult$pvalues  # observed p values
      simulatedCount = simulatedResult$simulated
    }
  }
  
  FDRResult = NA
  if (returnFDR) {
    # calculate the FDR for discrete counts
    observed = pvalues
    RBHResult = RBH(observed, nullP)
    
    FDRResult = cbind(observed, RBHResult)
    colnames(FDRResult) = c("pvalue", "RBH_point_estimate", "RBH_upper_limit")
    
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
    
    FDRResult = cbind(FDRResult, DBH, DBH.sd, ADBH.sd, data)
  }
  
  # get the mean of the ordered simulated pvalues
  nullPSorted = matrix(NA, dim(nullP)[1], dim(nullP)[2])
  for (j in 1:dim(nullP)[2]) {
    nullPSorted[, j] = sort(nullP[, j], decreasing = T)
  }
  nullPSortedMean = rowMeans(nullPSorted)
  
  lambdaEmp = NA
  intercept = NA
  if (!nolambda) {
    estimateResult = estimateLambdaEmp(pvalues, nullPSortedMean, quantile)
    lambdaEmp = estimateResult$lambdaEmp
    intercept = estimateResult$intercept
    if (is.na(lambdaEmp) || is.na(intercept)) { 
      print("cannot estimate lambda for qqplots, too few pvalues less than 1")
      plotNULL = F
    }
  }
  
  if (nolambda) {
    plotNULL = F
  }
  
  lambdaNull = NA
  if (plotNULL) {
    lambdaNull = numeric(nReplication)
    for (j in 1:nReplication) {
      nullPSortedMeanLOO = rowMeans(nullPSorted[, -j])
      lambdaNull[j] = estimateLambdaEmp(nullPSorted[, j], 
                                        nullPSortedMeanLOO, 0.95)[[1]]
    }
  }
  
  if (doPlot) {
    # truncate the p-value at 1e-30
    if (!is.null(limit)) {
      maxLogP = min(limit, maxLogP)
    }
    pvaluesPlot = pvalues
    pvaluesPlot[pvalues < 10^(-maxLogP)] = 10^(-maxLogP)
    observed = -log10(sort(pvaluesPlot, decreasing = T))
    expected = -log10(nullPSortedMean)
    
    if (is.null(limit)) {
      maxlog10 = ceiling(max(c(observed[length(observed)], 
                               expected[length(expected)])))
    } else {
      maxlog10 = limit
    }
    if (plotNULL) {
      par(fig=c(0,0.8,0,1))
    }
    plot(expected, observed,
         xlab = "Expected(-log10 P-value)", 
         ylab = "Observed(-log10 P-value)", type = "p", pch = 16, col = "blue",
         xlim = c(0, maxlog10), ylim = c(0, maxlog10), 
         cex.lab=cexRatio, cex.axis=cexRatio, cex.main=cexRatio, 
         cex.sub=cexRatio, ...
    )
    abline(a = 0, b = 1)
    if (!nolambda && !is.na(lambdaEmp) && !is.na(intercept)) {
      abline(a = intercept, b = lambdaEmp, lty = 3)
      #par(cex = cexRatio)
      legend("bottomright", 
             legend = substitute(paste(lambda[emp], '=', lambdaEmp, sep = ''),
             list(lambdaEmp = round(lambdaEmp, 2))),
             cex = cexRatio
      )
    }
    if (plotNULL) {
      par(fig=c(0.7,1,0,1), new=TRUE)
      boxplot(lambdaNull, xlab = expression(lambda[null]),
      #cex.lab=cexRatio, cex.axis=cexRatio, cex.main=cexRatio, cex.sub=cexRatio,
      cex = cexRatio)
    }
  } 
  
  if (!returnNullP) {
    nullP = NA
  }
  
  return(list(lambdaEmp = lambdaEmp, 
              lambdaNull = lambdaNull, 
              nullPSortedMean = nullPSortedMean,
              simulatedCount = simulatedCount, 
              nullP = nullP,
              pvalues = pvalues,
              FDRResult = FDRResult,
              indexKept = index,
              intercept = intercept))
  
}

#' QQ plot assuming uniform distribution of the P values
#'
#' It implements 3 different ways to estimate lambda
#' 1) twoPoints: based on the point with the observed pvalue of 1 with the 
#'     highest rank and the pvalue at the 95 percentile among the non-one 
#'     pvalues. Set twoPoints to TRUE to run this method
#' 2) 0.5-0.95: fit a line to the points within 0.5 and 0.95 quantiles.
#'     Set removeOne and twoPoints to FALSE to run this method
#' 3) removeOne: fit a line to the points starting with the point with 
#'   the observed pvalue of 1 and the highest rank and all points with 
#'   non-one pvalues until the point with the p-value of 0.95 quantile 
#'     Set removeOne to TRUE and twoPoints to FALSE to run this method. 
#' Reference: Gao et al, AJHG, 2018, Burden Testing of Rare Variants Identified 
#' through Exome Sequencing via Publicly Available Control Data. 

#' @param pvalues the pvalues used to estimate lambda and generate the QQ plot
#' @param doPlot whether to generate the QQ plot
#' @param upper the upper quantile of points used 
#' @param lower the lower quantile of points used
#' @param removeOne if true, remove all points with observed P value 1 except 
#' the point with the highest rank of expected values
#' @param twoPoints if ture, the lambda is estimated based on the point with 
#' the observed pvalue of 1 and the highest rank of expected p-values and 
#' the pvalue at the 95 percentile among the non-one 
#' @param cexRatio set the cex option in the plot
#' @param maxLogP the maximum of -log10(P) to be plotted, larger values will be
#' set to this value

#' @return the estimated lambda

#' @export
qqplotQT = function(pvalues, doPlot = T, upper = 0.95, lower = 0.5, 
                  removeOne = F, twoPoints = F, cexRatio = 1.5,
                  maxLogP = 30) {
  
  pvalues[pvalues < 10^(-maxLogP)] = 10^(-maxLogP)
  pvalues = pvalues[order(pvalues, decreasing = T)]
  expected = (length(pvalues):1) / (length(pvalues) + 1)
  
  log10P = -log10(pvalues)
  log10Expected = -log10(expected)
  
  maxlog10 = ceiling(max(c(log10P[length(log10P)], 
                           log10Expected[length(log10Expected)]), 
                         na.rm = T))
  
  regressionCoef = estimateLambdaQT(pvalues, upper = upper, lower = lower, 
                                    removeOne = removeOne, 
                                    twoPoints = twoPoints)
  a1 = regressionCoef[2]
  b1 = regressionCoef[1]
  
  if (doPlot) {
    # plot
    plot(log10Expected, log10P, xlab = "Expected(-log10 P-value)", 
         ylab = "Observed(-log10 P-value)", type = "p", pch = 16, col = "blue",
         xlim = c(0, maxlog10), ylim = c(0, maxlog10), 
    cex.lab=cexRatio, cex.axis=cexRatio, cex.main=cexRatio, cex.sub=cexRatio)
    abline(a = 0, b = 1)
    abline(a = a1, b = b1, lty = 3)
    #browser()
    #lambda = qchisq(median(pvalues), df = 1, lower.tail = F) / 0.455
    if (twoPoints) {
      legend("bottomright", 
             legend = substitute(paste(lambda["2points"], '=', b1, sep = ''),
                                 list(b1 = round(b1, 2))),
             cex = cexRatio
      )
    } else if (removeOne) {
      legend("bottomright", 
             legend = substitute(paste(lambda[remove1], '=', b1, sep = ''),
                                 list(b1 = round(b1, 2))),
             cex = cexRatio
      )
    } else {
      legend("bottomright", 
          legend = substitute(paste(lambda[lower-0.95], '=', b1, sep = ''),
                          list(b1 = round(b1, 2), lower = lower)),
          cex = cexRatio
      )
    }
  }
  
  return(b1)
}
