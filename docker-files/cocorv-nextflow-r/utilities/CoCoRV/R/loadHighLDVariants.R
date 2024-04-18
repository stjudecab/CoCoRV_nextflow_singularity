pairsToLists = function(variants) {
  # combine all pairwise LD variants into a list if connected components
  # If two variants in high LD, then they are connected
  
  edges = c(t(variants))
  g = make_undirected_graph(edges)
  
  highLDList = NULL
  connectedComponents = components(g)
  if (sum(connectedComponents$csize > 1)) {
    cluster = which(connectedComponents$csize > 1)
    for (i in 1:length(cluster)) {
      highLDList[[i]] = names(which(connectedComponents$membership == 
                                      cluster[i]))
    }
  }
  
  return(highLDList)
}

geneGroupedList = function(highLDVariant, ignoreEthnicityInLD) {
  highLDList = list()
  allGenes = unique(highLDVariant[, "gene"])
  for (i in 1:length(allGenes)) {
    index = (highLDVariant[, "gene"] == allGenes[i])
    LDPerGene = highLDVariant[index, , drop = F]
    if (ignoreEthnicityInLD) {
      highLDList[[i]] = pairsToLists(LDPerGene[, 1:2])
    } else {
      allEthnicities = unique(LDPerGene[, "ethnicity"])
      highLDList[[i]] = list()
      for (j in 1:length(allEthnicities)) {
        indexEthnicity = (LDPerGene[, "ethnicity"] == allEthnicities[j])
        highLDList[[i]][[j]] = pairsToLists(LDPerGene[indexEthnicity, 1:2])
        names(highLDList[[i]])[j] = allEthnicities[j]
      }
    }
    names(highLDList)[i] = allGenes[i]
  }
  
  return(highLDList)
}

#' load the variants in high LD. 

#' If the input file is not a rds file, the input file must have the 
#' following columns:
#' variant1  variant2  gene  ethnicity pvalue log10OR fdrEthnicity. If the 
#' input is a rds file, it will be loaded directly and other arguments will
#' not be applied

#' @param highLDVariantFile: a file containing the variants in high LD
#' @param ignoreEthnicityInLD: If true, ignore the ethnicity and merge all LDs in any 
#' @param ethnicity. Otherwise LDs are stored for each ethnicity
#' @param pThreshold: the p-value threshold for LD pairs 
#' @param ORThreshold: the odds ratio threshold for LD pairs

#' @return It returns a list of two compoents. The first one is based on 
#' the FDR threshold 0.05, and the second one is based on a more lenient 
#'  pvalue threshold specified by the argument pThreshold. Each component is also a list, where each element corresponds to a gene. 

#' @details For each gene, if ignoreEthnicityInLD is false, each gene has lists 
#' corresponding to ethnicities. For example, if ignoreEthnicityInLD is true, the 
#' highLDList for gene TAAR2 can be accessed as highLDList[["TAAR2"]]. If 
#' ignoreEthnicityInLD is false, it can be accessed as highLDList[["TAAR2"]][["afr"]] or 
#' highLDList[["TAAR2"]][["amr"]] 

#' @export
loadHighLDVariants = function(highLDVariantFile, ignoreEthnicityInLD = F,
                              pThreshold = 0.05,
                              ORThreshold = 1) { 
  
  if (regexpr(".rds$|.RDS$", highLDVariantFile) == -1) {
    highLDVariantFull = fread(highLDVariantFile, data.table = F)
    
    # stringent LD criterion based on fdr for use in controls to be conservative
    highLDFDR = highLDVariantFull[highLDVariantFull$fdrEthnicity < 0.05, , 
                                  drop = F]
    highLDListFDR = geneGroupedList(highLDFDR, ignoreEthnicityInLD)
    
    # relaxed LD criterion for use in cases with full genotypes 
    # for recessive models in order to be conservative
    # This also considers the fewer number of LD tests encountered when checking
    # variant pairs in cases
    highLDVariantRelaxed = highLDVariantFull[
                    highLDVariantFull$pvalue < pThreshold & 
                    highLDVariantFull$log10OR > log10(ORThreshold), , drop = F]
    highLDListRelaxed = geneGroupedList(highLDVariantRelaxed, 
                                        ignoreEthnicityInLD)
    return(list(highLDListFDR, highLDListRelaxed))
  } else {
    loadedLD = readRDS(highLDVariantFile)
    return(loadedLD)
  }
}



