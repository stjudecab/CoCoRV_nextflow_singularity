#' simulated null P values for both the FET and CMH through simulating
#' the actual 2*2*m tables under the null
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
#' @param alternative see alternative in fisher.test
#' @param nReplication the number of replications of empirical null P values 
#' @param returnSimulated whether to return simulated null count tables
#' @param returnNullP whether to return simulated null p values
#' @param verbose whether verbose or not
#' @param seed set the random seed

#' @return A list of the following components. 
#' simulated: the simulated count matrix under the null hypothesis
#' nullP: the matrix of null P values from the simulated counts
#' pvalues: observed p values

#' @export
simulateNullPvalues = function(data, alternative = "two.sided",
                    nReplication = 100, returnSimulated = F, returnNullP = T,
                    verbose = T, seed = NULL, ncore = 1) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  data = as.matrix(data)
  
  stopifnot(dim(data)[2] %% 4 == 0)
  nGroup = dim(data)[2] / 4
  
  # remove NA rows and rows with 0 mutations
  data = data[rowSums(is.na(data)) == 0 & 
            rowSums(data[, seq(from = 1, by = 2, to = dim(data)[2])]) > 0 & 
            rowSums(data[, seq(from = 2, by = 2, to = dim(data)[2])]) > 0, ]
  
  nGene = dim(data)[1]
  mM = matrix(NA, nGene, nGroup)
  nM = matrix(NA, nGene, nGroup)
  kM = matrix(NA, nGene, nGroup)
  
  for (i in 1:nGroup) {
    mM[, i] = data[, (i - 1) * 4 + 1] + data[, (i - 1) * 4 + 2]
    nM[, i] = data[, (i - 1) * 4 + 3] + data[, (i - 1) * 4 + 4]
    kM[, i] = data[, (i - 1) * 4 + 1] + data[, (i - 1) * 4 + 3]
  }
  
  x = array(NA, dim = c(dim(data)[1], nGroup, nReplication))
  
  # calculate the pvalues using exact tests
  pvalues = numeric(dim(data)[1])
  pvalues[] = NA
  for ( i in 1:dim(data)[1]) {
    if (nGroup == 1) {
      table = matrix(c(data[i, ]), nrow = 2, byrow = T)
    } else {
      # column first to fill the arrays, this should have the same results
      # as filling the row first
      table = array(data[i, ], dim = c(2, 2, nGroup))
      # make sure any strata with sample size > 1
      index = which(apply(table, 3L, sum) > 1)
      stopifnot(length(index) > 0)
      table = table[, , index, drop = F]
    }
    tryCatch({
      if (nGroup == 1) {
        pvalues[i] = fisher.test(table, alternative = alternative)$p.value
      } else {
        pvalues[i] = mantelhaen.test(table, exact = T, 
                                     alternative = alternative)$p.value
      }
    }, error = function(e) {
      print(e)
      print(table)
    })
    # if (is.na(pvalues[i])) {
    #   print("Use asympototic chi square distributions")
    #   if (nGroup == 1) {
    #     if (alternative == "two.sided") {
    #       pvalues[i] = chisq.test(table)$p.value
    #     } else {
    #       print(table)
    #       stop("Fisher's exact test failed")
    #     }
    #   } else {
    #     pvalues[i] = mantelhaen.test(table, exact = F,
    #                                  alternative = alternative)$p.value
    #   }
    # }
  }
  
  simulated = NA
  if (returnSimulated || returnNullP) {
    # simulate counts based on the null distribution
    simulated = list()
    for (i in 1:nGene) {
      for (j in 1:nGroup) {
        x[i, j, ] = rhyper(nReplication, mM[i, j], nM[i, j], kM[i, j])
      }
    }
    for (l in 1:nReplication) {
      simulatedMatrix = matrix(NA, nGene, nGroup * 4)
      for (j in 1:nGroup) {
        simulatedMatrix[, ((j - 1) * 4 + 1) : (j * 4)] = 
          cbind(x[, j, l], mM[, j] - x[, j, l], 
                kM[, j] - x[, j, l], nM[, j] - (kM[, j] - x[, j, l]))
      }
      simulated[[l]] = simulatedMatrix
    }
  }
  
  nullP = NA
  if (returnNullP) { # generate the null P values
    if (verbose) {
      cat(paste("Simulate empirical null P values:", nReplication, 
                "replicates\n"))
    }
    BPPARAM = initBPPARAM(ncore)
    nullP = bplapply(1:nReplication, calculateNullP, 
                  verbose = verbose, nReplication = nReplication, 
                  nGene = nGene, simulated = simulated, nGroup = nGroup,
                  alternative = alternative, BPPARAM = BPPARAM)
    nullP = matrix(unlist(nullP), nrow = nGene)
  }
  
  return(list(simulated = simulated, 
              nullP = nullP,
              pvalues = pvalues))
}

initBPPARAM = function(ncore) {
  if (ncore >= 1) {
    if (.Platform$OS.type == "windows") {
      parameter <- SnowParam(workers = ncore)
    } else {
      parameter <- MulticoreParam(workers = ncore)
    }
  } else {
    parameter = NULL
  }
  
  return(parameter)
}

calculateNullP = function(l, verbose, nReplication, nGene, simulated, 
                          nGroup, alternative) {
  if (verbose) {
    cat(".")
  }
  nullPVector = numeric(nGene)
  for (i in 1:nGene) {
    if (nGroup == 1) {
      table = matrix(simulated[[l]][i, ], nrow = 2, byrow = T)
    } else {
      table = array(simulated[[l]][i, ], dim = c(2, 2, nGroup))
      # make sure any strata with sample size > 1
      index = which(apply(table, 3L, sum) > 1)
      stopifnot(length(index) > 0)
      table = table[, , index, drop = F]
    }
    tryCatch({
      if (nGroup == 1) {
        nullPVector[i] = fisher.test(table, 
                                  alternative = alternative)$p.value
      } else {
        nullPVector[i] = mantelhaen.test(table, exact = T, 
                                      alternative = alternative)$p.value
      }
    },
    error = function(e) {
      print(e)
      print(table)
      stop()
    })
  }
  return(nullPVector)
}

#' simulated null P values for Fisher's exact test using the CDF of the 
#' p-values under the null directly. This is much faster than simulating the 
#' 2x2 tables as a intermediate step

#' @param data a  matrix or data frame with the number of columns being 4:
#' The columns should be arranged as follows
#'    counts with rare mutations in case, 
#'    counts without rare mutations in case, 
#'    counts with rare mutations in controls
#'    counts without rare mutations in controls
#' @param alternative see alternative in fisher.test
#' @param nReplication the number of replications of empirical null P values 
#' @param seed set the random seed

#' @return A list of the following components. 
#' pvalues: observed p values
#' nullP: the matrix of null P values from the simulated counts

#' @export
simulateNullPvaluesFET = function(data, alternative = "two.sided",
                                  nReplication = 1000,
                                  seed = NULL) {
  # simulate P values using the fast cdf based p-value sampling
  
  if (dim(data)[2] != 4) {
    print("The input must have four columns representing a 2x2 table")
    stop()
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  nGene = dim(data)[1]
  pvalueAndCDF = fisher.pvalues.support(counts = data, input = "noassoc", 
                                        alternative = alternative)
  observed = pvalueAndCDF$raw
  cdfList = pvalueAndCDF$support
  # simulate p-values directly from the cdfList
  simulated = matrix(NA, nGene, nReplication) 
  for (j in 1:nGene) {
    pmf = diff(c(0, cdfList[[j]])) 
    simulated[j, ] = sample(cdfList[[j]], size = nReplication, 
                            replace = T, prob = pmf)
  }
  
  return(list(pvalues = observed, nullP = simulated))
}

#' simulated null P values for CMH exact test using the CDF of the 
#' p-values under the null directly. This is much faster than simulating the 
#' 2x2xK tables as a intermediate step

#' @param data a  matrix or data frame with the number of columns being 4xK:
#' The columns of each 2x2 contingency table should be arranged as follows
#'    counts with rare mutations in case, 
#'    counts without rare mutations in case, 
#'    counts with rare mutations in controls
#'    counts without rare mutations in controls
#' @param alternative see alternative in fisher.test
#' @param nReplication the number of replications of empirical null P values 
#' @param seed set the random seed

#' @return A list of the following components. 
#' pvalues: observed p values
#' nullP: the matrix of null P values from the simulated counts

#' @export
simulateNullPvaluesCMH = function(data, alternative = "two.sided",
                                  nReplication = 1000,
                                  seed = NULL) {
  # simulate P values using the fast cdf based p-value sampling
  
  if (dim(data)[2] %% 4 != 0 || dim(data)[2] == 0) {
    print("The input must have 4xK columns representing 2x2xK tables")
    stop()
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }

  data = as.matrix(data)
  
  nGene = dim(data)[1]
  pvalueAndCDF = cmh.pvalues.support(counts = data,  
                                     alternative = alternative)
  observed = pvalueAndCDF$raw
  cdfList = pvalueAndCDF$support
  # simulate p-values directly from the cdfList
  simulated = matrix(NA, nGene, nReplication) 
  for (j in 1:nGene) {
    pmf = diff(c(0, cdfList[[j]])) 
    simulated[j, ] = sample(cdfList[[j]], size = nReplication, 
                            replace = T, prob = pmf)
  }
  
  return(list(pvalues = observed, nullP = simulated))
}