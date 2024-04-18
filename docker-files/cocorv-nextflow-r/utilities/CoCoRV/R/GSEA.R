#' Perform Gene set enrichment analysis using simulated null P values 
#' corresponding to phenotype permutation

#' @param data a  matrix or data frame with the number of columns being 4*k:
#'  If k is 1, it uses the Fisher exact test;
#'  if k > 1, it uses the CMH exact test. The columns should be arranged as
#'  follows for each stratified group
#'    counts with rare mutations in case, 
#'    counts without rare mutations in case, 
#'    counts with rare mutations in controls
#'    counts without rare mutations in controls
#' @param pathway a list of vectors with gene IDs for each pathway
#' @param alternative "two.sided", "greater", or "less". The later GSEA will 
#'     match the direction to perform standard GSEA, or one-sided GSEA.
#' @param p specify the power of the absolute value of the score, default 1.
#' @param nReplication number of actual replications of the null P values, 
#'    these will be used in later sampling to form permuted values in GSEA 
#' @param maxPermutation number of maximal permutations performed if the 
#'   maximal extreme value is not reached.
#' @param maxExtreme number of the maximal extreme values from the permutation 
#'    before stop, used to control the adaptive permutation

#' @return A data frame of GSEA results

#' @export
GSEA = function(data, pathway, alternative = "two.sided", p = 1, 
                minPathwaySize = 2,
                maxPathwaySize = Inf,
                nReplication = 1000, 
                maxPermutation = 1e6, 
                maxExtreme = 36, 
                ncore = 1) {
  
  # remove rows with 0 sums or NA
  indexRemove = rowSums(data[, seq(from = 1, by = 2, to = dim(data)[2]), 
                             drop = F]) == 0
  indexRemove[is.na(indexRemove)] = T
  
  if (sum(indexRemove) > 0) {
    data = data[!indexRemove, , drop = F]
    stopifnot(dim(data)[1] > 0)
  }
  
  # convert gene IDs in each path to gene index
  geneNames = rownames(data)
  pathwayIndex = list()
  k = 1
  for (i in 1:length(pathway)) {
    index = match(pathway[[i]], geneNames)
    index = index[!is.na(index)]
    if (length(index) >= minPathwaySize && length(index) <= maxPathwaySize) {
      pathwayIndex[[k]] = index
      names(pathwayIndex)[k] = names(pathway)[i]
      k = k + 1
    }
  }
  stopifnot(length(pathwayIndex) > 0)
  # randomize pathwayIndex for later parallel running
  pathwayIndex = pathwayIndex[sample(1:length(pathwayIndex))]
  
  tic = proc.time()
  
  # simulated = simulateNullPvalues(data, alternative = alternative,
  #                                       nReplication = nReplication, 
  #                                       returnSimulated = F, 
  #                                       returnNullP = T, ncore = ncore)
  
  if (dim(data)[2] > 4) {  ## CMH 
    # generate simualted p-values using CoCoRV
    # simuatedResult = simulateNullPvalues(data, nReplication = nReplication,
    #                                      alternative = alternative, 
    #                                      returnSimulated = F, 
    #                                      returnNullP = T,
    #                                      ncore = ncore)
    # faster sampling of p-values under null
    simuatedResult = simulateNullPvaluesCMH(data, nReplication = nReplication,
                                            alternative = alternative)
    observed = simuatedResult$pvalues
    simulated = simuatedResult$nullP
  } else {  ## FET 
    # generate the cdf list and sample the p-values directly
    simuatedResult = simulateNullPvaluesFET(data, nReplication = nReplication,
                                            alternative = alternative)
    observed = simuatedResult$pvalues
    simulated = simuatedResult$nullP
  }
  
  toc = proc.time()
  
  cat("Time used to simulate null P values\n")
  print(toc - tic)
  
  rm(simuatedResult)
  gc()
  
  scores = -log10(simulated)
  scoreObserved = -log10(observed)
  nGenes = dim(scores)[1]
  nReplicate = dim(scores)[2]
  
  scoresNull = c(scores)
  rm(scores)
  gc()
  
  # match the naming in fgsea
  if (alternative == "two.sided") {
    scoreType = "std"
  } else if (alternative == "greater") {
    scoreType = "pos"
  } else if (alternative == "less") {
    scoreType = "neg"
  }
  
  # find leading edges
  indexOrdered = order(scoreObserved, decreasing = T)
  scoreObservedSorted = scoreObserved[indexOrdered]
  leadingEdge = list()
  size = numeric(length(pathwayIndex))
  for (i in 1:length(pathwayIndex)) {
    inSet = logical(nGenes)
    inSet[ pathwayIndex[[i]] ] = T
    inSet = inSet[indexOrdered]
    ESResult = calculateESCppFast(scoreObservedSorted, inSet, scoreType, p, T)
    leadingEdge[[i]] = (geneNames[indexOrdered])[ ESResult[-1] + 1 ]
    size[i] = length(pathwayIndex[[i]])
  }
  
  # calculate p values
  if (ncore == 1) {
    pathwayResult = GSEACpp(pathway = pathwayIndex, scoresNull = scoresNull,
                            nGenes, nReplicate,
                            scoreObserved, scoreType,
                            maxExtreme, maxPermutation, p)
  } else if (ncore > 1) {
    BPPARAM = initBPPARAM(ncore)
    # generate a start-end list
    chunkSize = min(500, ceiling(length(pathwayIndex) / ncore))
    d = 1:length(pathwayIndex)
    indexList = split(d, ceiling(seq_along(d)/chunkSize))
    pathwayResultList = bplapply(indexList, calculatePvalues, 
                     pathwayIndex, scoresNull, 
                     nGenes, nReplicate,
                     scoreObserved, scoreType, 
                     maxExtreme, maxPermutation, p, 
                     BPPARAM = BPPARAM)
    # combine results
    for (i in 1:length(pathwayResultList)) {
      if (i == 1) {
        pathwayResult = pathwayResultList[[i]]
      } else {
        pathwayResult = rbind(pathwayResult, pathwayResultList[[i]])
      }
    }
  } else {
    stop(paste("wrong number of cores:", ncore))
  }
  
  colnames(pathwayResult) = c("ES", "pathwayPvalue", "nExtreme", "nTotal")
  
  # calculate FDR
  FDR = p.adjust(pathwayResult[, 2], method = "BH")
  
  pathwayID = names(pathwayIndex)
  pathwayResult = data.table(pathwayID, pathwayResult[, 1:2], FDR, 
                             pathwayResult[, 3:4], size, leadingEdge)

  return(pathwayResult)
}

calculatePvalues = function(index, pathwayIndex, scoresNull, 
                            nGenes, nReplicate,
                            scoreObserved, scoreType, 
                            maxExtreme, maxPermutation, p) {
  pathwayResult = GSEACpp(pathway = pathwayIndex[index], 
                          scoresNull = scoresNull, 
                          nGenes, nReplicate,
                          scoreObserved, scoreType, 
                          maxExtreme, maxPermutation, p)
  return(pathwayResult)
}