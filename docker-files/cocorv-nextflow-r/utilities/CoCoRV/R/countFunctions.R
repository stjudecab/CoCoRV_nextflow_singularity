calculateProb = function(p) {
  # calculate the probability of no alternate alleles, 1 alternate alleles, 
  # 2 alternate alleles and two heterogyzous genotypes
  
  # input
  # p: either a 1*nSNP matrix of AFs or a 3 * nSNP genotype frequencies
  # When p is a 3*nSNP matrix
  #   each row corresponds to genotype RR, RA, AA
  #   each column is a variant
  
  stopifnot(length(p) >= 1 || dim(p)[1] == 3 && dim(p)[2] >= 1)
  
  if (dim(p)[1] == 1) {
    geno = rbind((1 - p)^2, 2 * (1 - p) * p, p^2)
  } else {
    geno = p
  }
  nSNP = dim(geno)[2]
  
  if (nSNP == 1) {
    p0 = geno[1, 1]
    p1 = 1 - geno[1, 1] 
    p2 = geno[3, 1] 
    pTwoHets = 0
    p2h = p2
    pCHets = 0
    # p2Appr = p2
    # pTwoHetsAppr = pTwoHets
  } else {
    p0 = exp(sum(log(geno[1, ])))
    p1 = 0
    for (j in 1:nSNP) {
      genoj = geno[, j]
      genoi = geno[, -j, drop = F]
      p1 = p1 + exp(sum(log(genoi[1, ])) + log(1 - genoj[1]))
    }
    
    # exactly one homozyous alternate
    pOneHomAlt = 0
    for (j in 1:nSNP) {
      genoj = geno[, j]
      genoi = geno[, -j, drop = F]
      pOneHomAlt = pOneHomAlt + exp(log(genoj[3]) + sum(log(genoi[1, ])))
    }
    
    # two heterozygous 
    pTwoHets = 0
    if (nSNP == 2) {
      pTwoHets = geno[2, 1] * geno[2, 2]
    } else {
      for (i in 1:(nSNP - 1)) {
        for (j in (i + 1):nSNP) {
            genoi = geno[, i]
            genoj = geno[, j]
            genok = geno[, -c(i, j), drop = F]
            pTwoHets = pTwoHets + exp(log(1 - genoj[1]) + log(1 - genoi[1]) +
                                        sum(log(genok[1, ])))
            # <debug>  seems even smaller rmse based on simulation
            #pTwoHets = pTwoHets + exp(log(genoj[2]) + log(genoi[2]))
            # </debug>
        }
      }
    }
    
    p2 = pOneHomAlt + pTwoHets
    
    pCHets = pTwoHets / 2
    p2h = pOneHomAlt + pCHets
    
    # p2Appr = 1 - p0 - p1
    # pTwoHetsAppr = p2 - sum(p^2)
  }
  
  #return(cbind(p0, p1, p2, pTwoHets, p2Appr, pTwoHetsAppr))
  
  out = c(p0, p1, p2, pTwoHets, p2h, pCHets)
  names(out) = c("ZeroAlt", "OneAlt", "TwoAlt", "TwoHet", "TwoAltHap", "CHets")
  return(out)
}


#' estimate the counts of samples of different models based on AFs of variants
#' 
#' This function estimates the count of samples under the dominant model, the
#' recessive model and the double heterozygous model
#' 
#' @param nSample The total number of samples
#' @param AF The vector of alternate allele frequencies of each variant
#' @param excludeHighLD If true, exclude variants in high LD except the first
#'   variant with the largest AF when calculating the counts, otherwise 
#'   exclude all variants in high LD when calculating the counts and then 
#'   later add the counts due to the first variant of each LD group. 
#'   The default is T, so the variants with high LD are not considered in the 
#'   recessive or 2hets model.
#' @param highLDList A list of variants with high LD. Default NULL
#' @param countHomAlt Either NULL or a vector of the counts of the alternate
#'        homozygous genotypes. If not null, the estimation is based on 
#'        genotypes instead of assuming Hardy Weinberg equilibrium
#' @param twoHaps If true, two variants must be on the two different haplotypes
#'        to be considered as double heterozygous. Default F
#' 
#' @return A vector of three integer numbers correspond to the number of 
#' qualified samples under the dominant model, the recessive and the double
#' heterozygous model
#' 
#' @export
estimateCounts = function(nSample, AF, 
                          excludeHighLD = T,
                          highLDList = list(NULL),
                          countHomAlt = NULL,
                          twoHaps = F) {
  if (!is.null(countHomAlt)) {
    # update AF to genotype AF
    nHomAlt = countHomAlt
    nHet = nSample * 2 * AF - nHomAlt * 2
    nHet[nHet < 0] = 0 # avoid numerical error 
    pHomAlt = countHomAlt / nSample
    pHet = nHet / nSample
    pHomRef = 1 - pHomAlt - pHet
    inputAF = rbind(pHomRef, pHet, pHomAlt) 
  } else {
    inputAF = matrix(AF, nrow = 1)
  }
  
  highLDList2 = list(NULL)
  if (!is.null(highLDList[[1]])) {
    # keep only those in the AF
    highLDList3 = highLDList[[1]]
    j = 1
    for (i in 1:length(highLDList3)) {
      index = (highLDList3[[i]] %in% names(AF))
      if (sum(index) > 0) {
        highLDList2[[j]] = highLDList3[[i]][index]
        j = j + 1
      } 
    }
  }
  
  if (!is.null(highLDList2[[1]])) {
    firstVariantSet = sapply(highLDList2, 
                             FUN = function(x) {return(x[1])})
    NotFirstVariantSet = sapply(highLDList2, 
                                FUN = function(x) {return(x[-1])})
    
    # remove all variants with high LD except the first one 
    # (the one with the largest AF)
    used = !(names(AF) %in% NotFirstVariantSet)
    AFRecessive = inputAF[, used, drop = F]
    probs = calculateProbCpp(AFRecessive)
    dominantEstimate = round(nSample * (1 - probs[1]))
    
    if (twoHaps) {
      recessiveEstimate = round(nSample * probs[5])
      twoHetsEstimate = round(nSample * probs[6])
    } else {
      if (!excludeHighLD) {
        # remove all variants with high LD
        # and later add the counts of first variants with high LD as recessive 
        # and two hets
        used = !(names(AF) %in% unlist(highLDList2))
        AFRecessive = inputAF[, used, drop = F] 
        probs = calculateProbCpp(AFRecessive)
      }
      if (dim(AFRecessive)[2] > 0) {
        if (!twoHaps) {
          recessiveEstimate = round(nSample * probs[3])
          twoHetsEstimate = round(nSample * probs[4])
        } else {
          recessiveEstimate = round(nSample * probs[5])
          twoHetsEstimate = round(nSample * probs[6])
        }
      } else {
        recessiveEstimate = 0
        twoHetsEstimate = 0
      }
      if (!excludeHighLD) {
        # add the counts of the first variant (largest AF) with high LD  
        firstIndex = (names(AF) %in% firstVariantSet)
        if (!is.null(countHomAlt)) {
          frequencyOfRecessive = c(colSums(inputAF[2:3, firstIndex]), 
                                   probs[3]) 
          frequencyOf2Het = c(colSums(inputAF[2:3, firstIndex]), probs[4])
          #addedCounts = round(sum(inputAF[2:3, firstIndex]) * nSample)
        } else {
          frequencyOfRecessive = c(colSums(inputAF[2:3, firstIndex]), probs[3]) 
          frequencyOf2Het = c(colSums(inputAF[2:3, firstIndex]), probs[4])
          #addedCounts = round(sum(1 - (1 - AF[firstIndex])^2) * nSample)
        }
        
        recessiveEstimate = round((1 - prod(1 - frequencyOfRecessive)) *
                                    nSample)
        twoHetsEstimate = round((1 - prod(1 - frequencyOf2Het)) *
                                  nSample)
        
        # recessiveEstimate = recessiveEstimate + addedCounts
        # twoHetsEstimate = twoHetsEstimate + addedCounts
      }
    }
  } else {
    probs = calculateProbCpp(inputAF)
    dominantEstimate = round(nSample * (1 - probs[1]))
    if (!twoHaps) {
      recessiveEstimate = round(nSample * probs[3])
      twoHetsEstimate = round(nSample * probs[4])
    } else {
      recessiveEstimate = round(nSample * probs[5])
      twoHetsEstimate = round(nSample * probs[6])
    }
  }
  
  out = cbind(dominantEstimate, recessiveEstimate, twoHetsEstimate)
  
  return(out)
}

extractCountsFromGenotypes = function(genotypeMatrix, excludeHighLD = F, 
                                      highLDList = list(NULL)) {
  dominant = sum(rowSums(genotypeMatrix, na.rm = T) > 0, 
                 na.rm = T)
  if (!excludeHighLD || is.null(highLDList[[1]])) {
    recessive = sum(rowSums(genotypeMatrix, na.rm = T) >= 2, 
                    na.rm = T)
    twoHets = sum(rowSums(genotypeMatrix >= 1, na.rm = T) >= 2, 
                  na.rm = T)
  } else {
    highLDList = highLDList[[1]]
    pureRecessiveIndex = which(rowSums(genotypeMatrix == 2, na.rm = T) > 0)
    pureRecessive = length(pureRecessiveIndex)
    recessive = pureRecessive
    twoHets = 0
    index = which(rowSums(genotypeMatrix > 0, na.rm = T) >= 2)
    variantID = colnames(genotypeMatrix)
    if (length(index) > 0) {
      for (i in 1:length(index)) {
        variants = which(genotypeMatrix[index[i], ] > 0)
        notInLD = F
        for (j in 1:(length(variants) - 1)) {
          for (k in (j + 1):length(variants)) {
            vj = variants[j]
            vk = variants[k]
            # check whether in high LD
            inLD = F 
            for (l in 1:length(highLDList)) {
              if (variantID[vj] %in% highLDList[[l]] && 
                  variantID[vk] %in% highLDList[[l]]) {
                inLD = T
                break
              }
            }
            if (!inLD) {
              notInLD = T
              break
            }
          }
          if (notInLD) {
            break
          }
        }
        if (notInLD) {
          # avoid recount of those with pure recessive patterns
          if (!(index[i] %in% pureRecessiveIndex)) {
            recessive = recessive + 1
          }
          twoHets = twoHets + 1
        }
      }
    }
  }
  
  return(c(dominant, recessive, twoHets))
}

extractCounts = function(gds, sampleID, geneList, 
                         caseGroupInfo = NULL, groupIDSpecified = NULL,
                         excludeHighLD = F, highLDList = list(NULL)) {
  ## extract counts of three modes: 
  ##  dominant: count of individuals with at least one alternate alleles
  ##  recessive: count of individuals with at least two alternate alleles
  ##  twoHets: count of individuals with at least two heterozygous variants
  ##  additive: the count of alternate alleles
  
  ## input
  ## gds: the gds handle
  ## sampleID: a vector of sample IDs
  ## geneList: A list of groups of variants. For each component, 
  ##   the first column is the gene ID,
  ##   the rest are variant IDs. the variant IDs could be chr-pos-ref-alt or 
  ##   chr:pos_ref/alt or chr:pos_ref_alt or anything with separators ":", "|",
  ##   or "/". 
  ## caseGroupInfo: a vector of case group IDs of the sampleID, used for 
  ##   stratified analysis
  ## groupIDSpecified: a vector of the group IDs specified to match those in 
  #    controls
  ## excludeHighLD: if T, exclude variants in high LD when counting for 
  ##   the recessive or twoHets models
  
  ## output
  ## a list of two components
  #   extracted: a matrix of extracted cocunts with the following columns
  #      g1_dominant g2_dominant ... g1_recessive ... g1_twoHets ...
  #   estimated: if caseGroupInfo is not NULL, these are esitmated counts based
  #     on the extracted AC per group
  
  nGene = length(geneList)
  gene = character(nGene)
  
  # make sure the groupIDs are unique
  stopifnot(length(groupIDSpecified) == length(unique(groupIDSpecified)))
  
  nGroup = NULL
  estimated = NULL
  if (is.null(groupIDSpecified)) {
    nGroup = 1
  } else {
    nGroup = length(groupIDSpecified)
    estimated = matrix(0, nGene, nGroup * 3)
  }
  extracted = matrix(0, nGene, nGroup * 3)
  
  variantList = list()
  for (i in 1:nGene) {
    gene[i] = geneList[[i]][1]
    variantList[[i]] = gsub(":|_|/", "-", geneList[[i]][-1])
  }
  # convert to a long format
  numberPerGene = sapply(variantList, length)
  geneLongFormat = rep(gene, numberPerGene)
  variantsLongFormat = unlist(variantList)
  
  # get the variant IDs
  seqSetFilter(gds, verbose = F)
  variantsInData = seqGetData(gds, "$chrom_pos_allele")
  variantsInData = gsub(":|_", "-", variantsInData)
  
  indexLongFormat = match(variantsLongFormat, variantsInData)
  # convert to a index list
  indexList = split(indexLongFormat, 
                    factor(geneLongFormat, levels = unique(geneLongFormat)))
  
  for (i in 1:nGene) {
    selected = indexList[[i]]
    selected = selected[!is.na(selected)]
    if (length(selected) > 0) {
      #print(gene[i])
      for (j in 1:nGroup) {
        if (!is.null(groupIDSpecified) && !is.null(caseGroupInfo)) {
          sampleIDPerGroup = sampleID[caseGroupInfo == groupIDSpecified[j]]
        } else {
          sampleIDPerGroup = sampleID
        }
        if (length(sampleIDPerGroup) > 0) {
          seqSetFilter(gds, sample.id = sampleIDPerGroup, 
                       variant.sel = selected, 
                       verbose = F)
          genotypeMatrix = seqGetData(gds, "$dosage_alt")
          colnames(genotypeMatrix) = variantsInData[selected]
          if (!is.null(genotypeMatrix)) {
            #print(dim(genotypeMatrix))
            if (excludeHighLD) {
              counts = extractCountsFromGenotypes(genotypeMatrix, 
                                       excludeHighLD = excludeHighLD,
                                       highLDList = highLDList[i, j])
            } else {
              counts = extractCountsFromGenotypes(genotypeMatrix, 
                                          excludeHighLD = F,
                                          highLDList = list(NULL))
            }
            extracted[i, c(j, nGroup + j, nGroup * 2 + j)] = counts
            
            if (!is.null(estimated)) {
              nSample = length(sampleIDPerGroup)
              AC = colSums(genotypeMatrix, na.rm = T)
              AF = (AC + 0.1) / (2 * nSample)
              dominantEst = sum(AC)
              recessiveAndTwoHetsEst = estimateCounts(nSample, AF)[2:3]
              estimated[i, c(j, nGroup + j, nGroup * 2 + j)] = c(dominantEst,
                                                recessiveAndTwoHetsEst)
            }
          }
        }
      }
    }
  }
  
  out = list(extracted = extracted, estimated = estimated)
}

verifyInControlHighLD = function(caseHighLDList, controlHighLDList) {
  # make sure the high LD pairs are also in controls
  # If any pair of a LD set in case is also in the control, the whole set from
  # the case will be kept, and therefore later be excluded. 
  # This is a conservative strategy to avoid false positives
  
  # input
  # caseHighLDList, controlHighLDList: both are list of variant sets
  
  verified = NULL
  
  if (is.null(controlHighLDList[[1]]) || is.null(caseHighLDList[[1]])) {
    return(list(verified)) 
  }

  k = 1
  for (i in 1:length(caseHighLDList)) {
    keep = F
    for (j in 1:length(controlHighLDList)) {
      result = intersect(caseHighLDList[[i]], controlHighLDList[[j]])
      if (length(result) >= 2) {
        keep = T
        break
      }
    }
    if (keep) {
      verified[[k]] = caseHighLDList[[i]]
      k = k + 1
    }
  }
  
  return(list(verified))
}

extractHighLDFromCases = function(gds, sampleID, geneList, 
                         caseGroupInfo = NULL, groupIDSpecified = NULL,
                         ORThreshold,
                         pThreshold,
                         ACThreshold = 0,
                         controlHighLD = NULL,
                         checkHighLDInControl = F,
                         gnomADVersion = "v2exome",
                         controlCountQC = NULL,
                         groupColumn = NULL,
                         controlFunction = NULL,
                         controlACANConfig = NULL,
                         pLDControl = NULL,
                         ORThresholdControl = NULL,
                         highLDListRelaxed = NULL) {
  ## extract high LD from case genotype data. If checkHighLDInControl is T,
  ## it will also check whether variants are in LD in the controls
  
  ## input
  ## gds: the gds handle
  ## sampleID: a vector of sample IDs
  ## geneList: A list of groups of variants. For each component, 
  ##   the first column is the gene ID,
  ##   the rest are variant IDs. the variant IDs could be chr-pos-ref-alt or 
  ##   chr:pos_ref/alt or chr:pos_ref_alt or anything with separators ":", "|",
  ##   or "/". 
  ## caseGroupInfo: a vector of case group IDs of the sampleID, used for 
  ##   stratified analysis
  ## groupIDSpecified: a vector of the group IDs specified to match those in 
  #    controls
  
  ## output
  ## a gene * group matrix, each element is a list of high LD variants
  
  nGene = length(geneList)
  genes = character(nGene)
  
  # make sure the groupIDs are unique
  stopifnot(length(groupIDSpecified) == length(unique(groupIDSpecified)))
  
  nGroup = NULL
  estimated = NULL
  if (is.null(groupIDSpecified)) {
    nGroup = 1
  } else {
    nGroup = length(groupIDSpecified)
    estimated = matrix(0, nGene, nGroup * 3)
  }
  extracted = matrix(0, nGene, nGroup * 3)
  
  variantList = list()
  for (i in 1:nGene) {
    genes[i] = geneList[[i]][1]
    variantList[[i]] = gsub(":|_|/", "-", geneList[[i]][-1])
  }
  # convert to a long format
  numberPerGene = sapply(variantList, length)
  geneLongFormat = rep(genes, numberPerGene)
  variantsLongFormat = unlist(variantList)
  
  # get the variant IDs
  seqSetFilter(gds, verbose = F)
  variantsInData = seqGetData(gds, "$chrom_pos_allele")
  variantsInData = gsub(":|_", "-", variantsInData)
  
  indexLongFormat = match(variantsLongFormat, variantsInData)
  # convert to a index list
  indexList = split(indexLongFormat, 
                    factor(geneLongFormat, levels = unique(geneLongFormat)))
  
  highLDListAll = matrix(rep(list(NULL), nGene * nGroup), nGene, nGroup)
  
  # for add LD in controls
  controlDataForTest = NULL
  ethnicitySet = ""
  
  for (i in 1:nGene) {
    selected = indexList[[i]]
    selected = selected[!is.na(selected)]
    if (length(selected) >= 2) {
      gene = genes[i]
      if (checkHighLDInControl && is.null(highLDListRelaxed)) {
        ## set up control data
        variantInControlSelected = (controlCountQC[, groupColumn] == gene & 
                                      controlFunction)
      }
      
      for (j in 1:nGroup) {
        if (!is.null(groupIDSpecified) && !is.null(caseGroupInfo)) {
          sampleIDPerGroup = sampleID[caseGroupInfo == groupIDSpecified[j]]
        } else {
          sampleIDPerGroup = sampleID
        }
        if (length(sampleIDPerGroup) > 0) {
          seqSetFilter(gds, sample.id = sampleIDPerGroup, 
                       variant.sel = selected, 
                       verbose = F)
          genotypeMatrix = seqGetData(gds, "$dosage_alt")
          colnames(genotypeMatrix) = variantsInData[selected]
          # keep variants with genotypes missingness < 0.5 and AF > 0 & AF < 1
          missingFreq = colMeans(is.na(genotypeMatrix))
          AF = colMeans(genotypeMatrix, na.rm = T) / 2
          used = (missingFreq < 0.5 & AF > 0 & AF < 1)
          genotypeMatrix = genotypeMatrix[, used, drop = F]
          if (!is.null(genotypeMatrix) && dim(genotypeMatrix)[1] > 1 &&
              dim(genotypeMatrix)[2] > 1) {
            #print(dim(genotypeMatrix))
            
            if (checkHighLDInControl) {
              if (is.null(highLDListRelaxed)) {
                controlACGroup = controlACANConfig[j, 1]  
                ethnicitySetAll = c("nfe", "afr", "amr", "eas", "sas", "fin")
                if (controlACGroup == "AC") {
                  ethnicitySet = ethnicitySetAll
                } else {
                  ethnicityIndex = grepl(sub(".*AC.*_", "", controlACGroup), 
                                         ethnicitySetAll)
                  if (sum(ethnicityIndex) == 0) {
                    stop(paste(
                      "cannot identify the ethnicity group from AC ID:",
                               controlACGroup))
                  }
                  ethnicitySet = ethnicitySetAll[ethnicityIndex]
                }
                features = c(1, grep("(AC|AN).*(nfe|amr|afr|eas|sas|fin)", 
                                     colnames(controlCountQC)))
                controlDataPerGene = controlCountQC[variantInControlSelected, 
                                                    features, 
                                                    drop = F]
                controlDataForTest = extractIndependentGnomADCounts(
                                                          controlDataPerGene,
                                                          ethnicitySet,
                                                          gnomADVersion)
              } else {
                controlDataForTest = NULL
              }
            }
            
            highLDList = extractHighLD(genotypeMatrix, 
                               ORThreshold,
                               pThreshold,
                               ACThreshold,
                               controlDataForTest = controlDataForTest,
                               pLDControl = pLDControl,
                               ORThresholdControl = ORThresholdControl,
                               highLDListRelaxed = highLDListRelaxed[i, j],
                               ethnicitySet = ethnicitySet)
            if (is.null(controlHighLD) || is.null(highLDList[[1]])) {
              highLDListAll[i, j] = highLDList
            } else {
              highLDListAll[i, j] = verifyInControlHighLD(highLDList[[1]], 
                                                  controlHighLD[i, j][[1]])
            }
          }
        }
      }
    }
  }
  
  return(highLDListAll)
}

extractHighLD = function(genotypeMatrix, ORThreshold = log10(1e3), 
                         pThreshold = 0.05, ACThreshold = 0,
                         controlDataForTest = NULL,
                         pLDControl = NULL,
                         ORThresholdControl = NULL,
                         highLDListRelaxed = NULL,
                         ethnicitySet = NULL) {
  ## extract variants of high LD from the genotype matrix
  
  ## input
  ## genotypeMatrix: the genotype matrix of size nSample * nVariant
  
  ## output
  ## a list of groups of variants with high LD
  
  stopifnot(dim(genotypeMatrix)[1] > 1 & dim(genotypeMatrix)[2] > 1)
  highLDList = NULL
  
  # check variant pairs that are in the same individual
  indWithMultipleVariants = which(rowSums(genotypeMatrix > 0, na.rm = T) >= 2)
  if (length(indWithMultipleVariants) == 0) {
    return(list(highLDList))
  }
  
  # order variants by MAF
  AF = colMeans(genotypeMatrix, na.rm = T) / 2
  MAF = AF
  MAF[AF > 0.5] = 1 - AF[AF > 0.5]
  genotypeMatrix = genotypeMatrix[, order(MAF)]
  
  nVariant = dim(genotypeMatrix)[2]
  LDMatrix = matrix(NA, nVariant, nVariant)
  colnames(LDMatrix) = colnames(genotypeMatrix)
  rownames(LDMatrix) = colnames(LDMatrix)
  for (i in 1:length(indWithMultipleVariants)) {
    variants = which(genotypeMatrix[indWithMultipleVariants[i], ] > 0)
    for (j in 1:(length(variants) - 1)) {
      for (k in (j + 1):length(variants)) {
        vj = variants[j]
        vk = variants[k]
        if (is.na(LDMatrix[vj, vk])) {
          if (pThreshold > 0) {
            testResult = LDTest(
              genotypeMatrix[, c(vj, vk)], 
              ORThreshold = ORThreshold,
              alternative = "greater",
              genotype = T,
              ACThreshold = ACThreshold)
          } else {
            testResult = 1
          }
          inLDInControl = F
          if (testResult[1] >= pThreshold) { 
            # check the control LD and add it if in LD
            vjID = colnames(genotypeMatrix)[vj]
            vkID = colnames(genotypeMatrix)[vk]
            if (!is.null(highLDListRelaxed)) {
              # directly look up the precomputed LD list
              for (l in 1:length(highLDListRelaxed)) {
                if (vjID %in% highLDListRelaxed[[l]] && 
                    vkID %in% highLDListRelaxed[[l]]) {
                  inLDInControl = T
                  break
                }
              }
            } else {
              if (!is.null(controlDataForTest)) {
                matchIndex = match(c(vjID, vkID), 
                                   colnames(controlDataForTest[[1]])) 
                if (!is.na(matchIndex[1]) && !is.na(matchIndex[2])) {
                  summaryControlCount = controlDataForTest[[1]][, matchIndex]
                  dataForEstimation = cbind(summaryControlCount, 
                                            controlDataForTest[[2]])
                  resultInControl = LDTest(
                                           dataForEstimation, 
                                           ORThresholdControl,
                                           ACThreshold = 0)
                  if (!is.na(resultInControl[1]) && 
                      resultInControl[1] < pLDControl) {
                    inLDInControl = T
                  }
                }
              }
            }
          }
          if (testResult[1] < pThreshold || inLDInControl) {
            LDMatrix[vj, vk] = T
            LDMatrix[vk, vj] = T
          } else {
            LDMatrix[vj, vk] = F
            LDMatrix[vk, vj] = F
          }
        } 
      }
    }
  }
  LDMatrix[is.na(LDMatrix)] = F
  
  # convert to a list
  g  = graph.adjacency(LDMatrix)
  connectedComponents = components(g)
  if (sum(connectedComponents$csize > 1)) {
    cluster = which(connectedComponents$csize > 1)
    for (i in 1:length(cluster)) {
      highLDList[[i]] = colnames(LDMatrix)[
           which(connectedComponents$membership == cluster[i])]
    }
  } 
  
  return(list(highLDList))
}

#' This function calculate AC AN from a gds File stratified in groups

#' @param gds a gds file handle
#' @param sampleID restrict to the samples specified 
#' @param bedGRange If set, it is a GRanges or GRangesList to specify a region
#' @param caseGroupInfo A vector of case group ID
#' @param groupIDSpecified A vector of group ID used
#' @param overlapType the overlap type between the bedGRange and the variants
#' in the gdsFile. See type in the findOverlaps function in GenomicRanges
#' @param nVariant number of variants
#' @param sexID the sex ID, 1 for male, 2 for female, used for counting the
#' total alleles in chromosome X and Y, default is NULL 
#' @param reference the reference build, used to identify the PAR region in 
#' chromosome X and Y, default is NULL


#' @return a list of two components. The first is the AC matrix of different
#' groups, the second is the AN matrix of different groups 

#' @export
extractACANCounts = function(gds, sampleID, bedGRange,
                            caseGroupInfo, groupIDSpecified, nVariant,
                            overlapType = "within", 
                            sexID = NULL, 
                            reference = NULL) {
  nGroup = length(groupIDSpecified)
  ACMatrix = matrix(0L, nVariant, nGroup)
  ANMatrix = matrix(0L, nVariant, nGroup)
  stopifnot(length(sampleID) > 0 && length(groupIDSpecified) > 0)
  for (i in 1:nGroup) {
    sampleIDInGroup = sampleID[caseGroupInfo == groupIDSpecified[i]]
    sexIDInGroup = NULL
    if (!is.null(sexID)) {
      sexIDInGroup = sexID[caseGroupInfo == groupIDSpecified[i]]
    }
    if (length(sampleIDInGroup) > 0) {
      
      # set filters
      seqSetFilter(gds, verbose = F)
      if (overlapType == "gds") {
        seqSetFilter(gds,  variant.sel = bedGRange, verbose = F)
        if (!is.null(sampleIDInGroup)) {
          seqSetFilter(gds,  sample.id = sampleIDInGroup, action = "intersect", 
                       verbose = F)
        }
      } else {
        if (!is.null(bedGRange)) {
          # get chr pos alleles
          chromosome = seqGetData(gds, "chromosome")
          position = seqGetData(gds, "position")
          alleles = seqGetData(gds, "allele")
          allelesMatrix = matrix(unlist(strsplit(alleles, ",")), ncol = 2, 
                                 byrow = T)
          allelesMatrix = as.matrix(allelesMatrix)
          lengthRef = nchar(allelesMatrix[, 1])
          lengthAlt = nchar(allelesMatrix[, 2])
          alleleMaxLength = lengthAlt
          indexRefLarger = (lengthRef > lengthAlt)
          alleleMaxLength[indexRefLarger] = lengthRef[indexRefLarger]
          
          variantGRange = GRanges(seqnames = chromosome,
                                  ranges = IRanges(start = position,
                                        end = position + alleleMaxLength - 1),
                              strand = Rle(strand("+"), length(chromosome)))
          findCovered = findOverlaps(variantGRange, bedGRange, 
                                     type = overlapType, select = "first")
          
          selectionIndex = !is.na(findCovered)
        } else {
          selectionIndex = NULL
        }
        seqSetFilter(gds, sample.id = sampleIDInGroup, 
                     variant.sel = selectionIndex,
                     verbose = F)
      }
      
      AC = seqApply(gds, "$dosage_alt", as.is="double", margin="by.variant",
                    FUN=function(x) { sum(x, na.rm=TRUE) })
      # check whether there are variants in X or Y non-PAR
      chromosome = seqGetData(gds, "chromosome")
      position = seqGetData(gds, "position")
      isXYNonPAR = XYNonPAR(chromosome, position, reference)
      AN = seqApply(gds, "$dosage_alt", as.is="double", margin="by.variant",
                    FUN=function(x) { sum(!is.na(x)) * 2 })
      if (sum(isXYNonPAR) > 0) {
        # for variants in XYNonPAR
        # I assume samples in each group have the same sex
        # for females, 
        # if chr == X, it is diploid, if chr == Y, it is 0
        # for males, it is haploid
        stopifnot(!is.null(sexIDInGroup))
        if (sexIDInGroup[1] == 1) { # male group
          AN[isXYNonPAR] = AN[isXYNonPAR] / 2
        } else { # female group
          AN[isXYNonPAR & grepl("Y", chromosome)] = 0
        } 
      }

      ACMatrix[, i] = AC
      ANMatrix[, i] = AN
    }
  }
  
  return(list(AC = ACMatrix, AN = ANMatrix))
}

#' This function extracts annotations from a gds File

#' @param gds a gds file handle
#' @param sampleID If set, then will restrict to the samples specified first
#' @param bedGRange If set, it is a GRanges or GRangesList to specify a region
#' @param updateACAN whether to calculate the AC AN instead of extracting from
#' the annotations
#' @param overlapType the overlap type between the bedGRange and the variants
#' in the gdsFile. See type in the findOverlaps function in GenomicRanges

#' @return a matrix with variantID as the first column and the rest are 
#' annotations in the gds file

#' @export
extractACANAnnotations = function(gds, sampleID = NULL, 
                                  bedGRange = NULL, 
                                  updateACAN = F,
                                  overlapType = "within",
                                  ACANIDs = NULL, 
                                  annotationSet = NULL) {
  # Note: those annotations with variable length will not be extracted. These
  #  annotations are not likely to be useful because of the variable length
  #  among variants
  
  seqSetFilter(gds, verbose = F)
  if (overlapType == "gds") {
    seqSetFilter(gds,  variant.sel = bedGRange, verbose = F)
    if (!is.null(sampleID)) {
      seqSetFilter(gds,  sample.id = sampleID, action = "intersect", 
                 verbose = F)
    }
  } else {
    if (!is.null(bedGRange)) {
      # get chr pos alleles
      chromosome = seqGetData(gds, "chromosome")
      position = seqGetData(gds, "position")
      # get the max length of alleles
      # alleleMaxLength <- seqApply(gds, "allele",
      #         FUN=function(x) max(nchar(unlist(strsplit(x,",")))),
      #         as.is="integer")
      
      alleles = seqGetData(gds, "allele")
      allelesMatrix = matrix(unlist(strsplit(alleles, ",")), ncol = 2, byrow = T)
      allelesMatrix = as.matrix(allelesMatrix)
      lengthRef = nchar(allelesMatrix[, 1])
      lengthAlt = nchar(allelesMatrix[, 2])
      alleleMaxLength = lengthAlt
      indexRefLarger = (lengthRef > lengthAlt)
      alleleMaxLength[indexRefLarger] = lengthRef[indexRefLarger]
      
      variantGRange = GRanges(seqnames = chromosome,
                              ranges = IRanges(start = position,
                                          end = position + alleleMaxLength - 1),
                              strand = Rle(strand("+"), length(chromosome)))
      findCovered = findOverlaps(variantGRange, bedGRange, 
                                 type = overlapType, select = "first")
      
      selectionIndex = !is.na(findCovered)
    } else {
      selectionIndex = NULL
    }
    seqSetFilter(gds, sample.id = sampleID, variant.sel = selectionIndex,
                   verbose = F)
  }

  # get the variant IDs
  variantID = seqGetData(gds, "$chrom_pos_allele")
  if (length(variantID) == 0) {
    return(NULL)
  }
  
  variantID = gsub(":|_", "-", variantID)
  FILTER = seqGetData(gds, "annotation/filter")
  
  # get annotations
  annotationInfo = seqSummary(gds, "annotation/info", verbose = F)
  
  if (updateACAN) {
    AC = seqApply(gds, "$dosage_alt", as.is="double", margin="by.variant",
                  FUN=function(x) { sum(x, na.rm=TRUE) })
    AN = seqApply(gds, "$dosage_alt", as.is="double", margin="by.variant",
                  FUN=function(x) { sum(!is.na(x)) * 2 })
    # extract all annotations except AC, AN
    annotationID = setdiff(annotationInfo$ID, c("AC", "AN"))
  } else {
    annotationID = annotationInfo$ID
  }
  
  # do not extract some unrelated annotations 
  removeIndex = which(grepl(
        "RAW_MQandDP|vep|hist|faf", 
                       annotationID))
  keepIndex = which(grepl("popmax", annotationID))
  removeIndexFinal = setdiff(removeIndex, keepIndex)
  if (length(removeIndexFinal) > 0) {
    annotationID = annotationID[-(removeIndexFinal)]
  }

  # if annotationSet is specified, load only those specified ones plus the 
  # ACANIDs if needed
  if (!is.null(annotationSet)) {
    if (!updateACAN) {  
    # this means full genotype data is available, so no need to include 
    # ACANIDs
      annotationSet = union(ACANIDs, annotationSet)
    } else {
      annotationSet = setdiff(annotationSet, c("AC", "AN"))
    }
    # check whether all specified annotations are in the data
    if (length(annotationSet) > 0) {
      insideData = (annotationSet %in% annotationID)
      if (!all(insideData)) {
        print("The following annotations are not found in the data")
        print(annotationSet[!insideData])
        stop()
      } else {
        annotationID = annotationSet
      }
    } else {
      annotationID = NULL  # this could happen if updateACAN is T
    }
  }
  
  if (length(annotationID) > 0) {
    annotationData = list()
    k = 1
    extractedIDs = NULL
    for (i in 1:length(annotationID)) {
      print(paste("extract annotation", i, ":", annotationID[i]))
      annotation = seqGetData(gds, paste0("annotation/info/", annotationID[i]))
      if (class(annotation)[1] == "SeqVarDataList") {
        # deal with SeqVarDataList
        # print(annotation[1])
        # browser()
        annotation = seqListVarData(annotation)
        annotation = sapply(annotation, function(x) {
                              if (length(x) > 1)
                                return(paste(x, collapse = ","))
                              else if (length(x) == 0) {
                                return(NA)
                              } else 
                                return(x)
                              })
      }
      annotationData[[k]] = annotation
      k = k + 1
      extractedIDs = c(extractedIDs, annotationID[i])
    }
    
    #browser()
    
    annotationData = do.call(cbind.data.frame, annotationData)
    names(annotationData) = extractedIDs
    
    if (updateACAN) {
      ACANAnnotation = cbind(variantID, FILTER, AC, AN, annotationData)
    } else {
      ACANAnnotation = cbind(variantID, FILTER, annotationData)
    }
  } else {
    if (updateACAN) {
      ACANAnnotation = data.frame(variantID, FILTER, AC, AN)
    } else {
      stop("no annotations availabe in the GDS file")
    }
  }
  
  return(ACANAnnotation)
}

generateCountTable = function(caseCount, nCase, controlCount, nControl, 
                              groupIDs = "") {
  # generate a count matrix where for each group, it is of 4 numbers: 
  #  caseMutation caseWOMutation controlMutation controlWOMutation
  # For all data, row represents genes, column represents groups
  
  # input
  # caseCount: a matrix of counts with variants of interest in cases 
  # nCase: a matrix of counts, the total sample size per group of cases
  # controlCount: a matrix of counts with variants of interest in controls
  # nControl: a matrix of counts, the total sample size per group of controls
  
  # output
  # a matrix of counts of 4 * nGroup, when nGroup = 1, can be used for Fisher
  # exact test, nGroup > 1, CMH exact test, when converint to a matrix or array
  
  nGene = dim(nCase)[1]
  nGroup = dim(nCase)[2]
  
  stopifnot(nGene > 0 && nGroup >= 0 && 
              identical(dim(caseCount), dim(controlCount)) && 
              identical(dim(nCase), dim(nControl)))
  
  countMatrix = matrix(NA, nGene, 4 * nGroup)
  columnNames = c("caseWMutation", 
                  "caseWOMutation",
                  "controlWMutation",
                  "controlWOMutation")

  countMatrix[, seq(from = 1, by = 4, to = 4 * nGroup)] = 
                                   caseCount
  countMatrix[, seq(from = 2, by = 4, to = 4 * nGroup)] = 
                            nCase - caseCount
  countMatrix[, seq(from = 3, by = 4, to = 4 * nGroup)] = 
                            controlCount
  countMatrix[, seq(from = 4, by = 4, to = 4 * nGroup)] = 
                           nControl - controlCount
  colnames(countMatrix) = rep("", dim(countMatrix)[2])
  for (i in 1:nGroup) {
    colnames(countMatrix)[((i - 1) * 4 + 1) : (i * 4)] = 
            paste0(columnNames, groupIDs[i])
  }
  
  return(countMatrix)
}

#' Perform contingency table based tests
#' 
#' If the number of columns is 4, Fisher exact test is used
#' If the number of columns is a multiple of 4 (> 4), CMH exact test is used

#' @param countMatrix a matrix of 4 * nGroup, ordered as
#'  caseWithMutation caseWOMutation controlWithMutation controlWOMutation
#'  and then by groups
#' @param test name of the test. "FET" means Fisher exact test, "CMH" means CMH
#   exact test

#' @return a matrix of four columns: p-value, odds ratio, 
#'   95% confidence interval lower bound,
#'   95% confidence interval upper bound

#' @export
contingencyTableTest = function(countMatrix, test = "FET") {
  nGene = dim(countMatrix)[1]
  stopifnot(nGene > 0 && dim(countMatrix)[2] %% 4 == 0)
  nGroup = dim(countMatrix)[2] / 4
  
  result = matrix(NA, nGene, 4)
  for (i in 1:nGene) {
    if (test == "FET" && dim(countMatrix)[2] == 4) {
      tryCatch({
        testResult = fisher.test(matrix(countMatrix[i, ], 2, 2))
        result[i, 1] = testResult$p.value
        result[i, 2] = testResult$estimate
        result[i, 3:4] = testResult$conf.int
      }, error = function(e) {
        print(countMatrix[i, ])
        print(e)
      })
    } else if (test == "CMH" && nGroup > 1) {
      data = array(countMatrix[i, ], dim = c(2, 2, nGroup))
      tryCatch({
        testResult = mantelhaen.test(data, exact = T)
        result[i, 1] = testResult$p.value
        result[i, 2] = testResult$estimate
        result[i, 3:4] = testResult$conf.int
      }, error = function(e) {
        print(data)
        print(e)
      })
    } else {
      stop("no test can be applied")
    }
  }
  
  colnames(result) = c("P", "OR", "CI_Lower", "CI_Upper")
  return(result)
}

#' extract the independent counts from gnomAD counts

#' @param controlDataPerGene a matrix, where the first column is the variant ID, 
#' and the rest are gnomAD annotated counts.
#' @param ethnicitySet which ethnicities to use in the extraction. For example, 
#' the 6 major ethnicities: c("nfe", "afr", "amr", "eas", "sas", "fin")
#' @param gnomADVersion gnomAD major version (v2exome, v2genome or v3genome), 
#' default is v2exome. 

#' @return a list of two components. The first component is a matrix, where
#' each column is a vairant, and rows are independent set of summary counts. 
#' The second component is a vector of averaged ANs across the variants, 
#' representing an estimate of the total number of haplotyppes pooled.

#' @export
extractIndependentGnomADCounts = function(controlDataPerGene,
                                          ethnicitySet, 
                                          gnomADVersion = "v2exome") {
  stopifnot(gnomADVersion %in% c("v2exome", "v2genome", "v3genome"))
  data = NULL
  
  if (gnomADVersion %in%  c("v2exome", "v2genome")) {
    for (j in 1:length(ethnicitySet)) {
      ethnicity = ethnicitySet[j] 
      if (ethnicity != "") {
        ACString = paste0("AC_", ethnicity)
        ANString = paste0("AN_", ethnicity)
      } else {
        ACString = "AC"
        ANString = "AN"
      }
      features = grep(paste0("^", ACString, "|_", ACString, "|", 
                             "^", ANString, "|_", ANString), 
                      colnames(controlDataPerGene), value = "T")
      if (is.null(data)) {
        data = controlDataPerGene[, features, drop = F]
        rownames(data) = controlDataPerGene[, 1]
      } else {
        data = cbind(data, controlDataPerGene[, features, drop = F])
      }
      
      # use uncorrelated counts
      gender = c("male", "female")

      for (i in 1:length(gender)) {
        if (gnomADVersion == "v2exome") {
          # cancer 
          feature = paste0("cancer_", ACString, "_", gender[i])
          data[, feature] = data[, paste0(ACString, "_", gender[i])] - 
            data[, paste0("non_", feature)]
          
          feature = paste0("cancer_", ANString, "_", gender[i])
          data[, feature] = data[, paste0(ANString, "_", gender[i])] - 
            data[, paste0("non_", feature)]
          
          # non_cancer_non_control
          feature = paste0("non_cancer_non_control_", ACString, "_", gender[i])
          data[, feature] = 
            data[, paste0("non_cancer_", ACString, "_", gender[i])] - 
            data[, paste0("controls_", ACString, "_", gender[i])]
          
          feature = paste0("non_cancer_non_control_", ANString, "_", gender[i])
          data[, feature] = 
            data[, paste0("non_cancer_", ANString, "_", gender[i])] - 
            data[, paste0("controls_", ANString, "_", gender[i])]
        } else { # v2gnome
          feature = paste0("non_control_", ACString, "_", gender[i])
          data[, feature] = 
            data[, paste0(ACString, "_", gender[i])] - 
            data[, paste0("controls_", ACString, "_", gender[i])]
          
          feature = paste0("non_control_", ANString, "_", gender[i])
          data[, feature] = 
            data[, paste0(ANString, "_", gender[i])] - 
            data[, paste0("controls_", ANString, "_", gender[i])]
        }
      }
    }
  
    # independent counts
    if (gnomADVersion == "v2exome") {
      independentFeatures = grep(
    "^cancer_AC_.*male|^non_cancer_non_control_AC_.*male|^controls_AC_.*male", 
      colnames(data), value = T)
    } else {
      independentFeatures = grep("^non_control_AC_.*male|^controls_AC_.*male", 
        colnames(data), value = T)
    }
    independentFeatures = sort(independentFeatures)
    independentFeaturesAN = sub("_AC_", "_AN_", independentFeatures)
  } else {
    for (j in 1:length(ethnicitySet)) {
      ethnicity = ethnicitySet[j] 
      if (ethnicity != "") {
        ethnicityString = paste0("_", ethnicity)
      } else {
        ethnicityString = ""
      }
      features = grep(paste0("^AC.*", ethnicityString, 
                             "|^AN.*", ethnicityString), 
                      colnames(controlDataPerGene), value = "T")
      if (is.null(data)) {
        data = controlDataPerGene[, features, drop = F]
        rownames(data) = controlDataPerGene[, 1]
      } else {
        data = cbind(data, controlDataPerGene[, features, drop = F])
      }
      
      # use uncorrelated counts
      gender = c("XY", "XX")
     
      for (i in 1:length(gender)) {
        # cancer 
        feature = paste0("AC_cancer", ethnicityString, "_", gender[i])
        data[, feature] = data[, 
                               paste0("AC", ethnicityString, "_", gender[i])] - 
          data[, paste0("AC_non_cancer", ethnicityString, "_", gender[i])]
        
        feature = paste0("AN_cancer", ethnicityString, "_", gender[i])
        data[, feature] = data[, 
                               paste0("AN", ethnicityString, "_", gender[i])] - 
          data[, paste0("AN_non_cancer", ethnicityString, "_", gender[i])]
        
        # non_cancer_non_control
        feature = paste0("AC_non_cancer_non_control", ethnicityString, "_", 
                         gender[i])
        data[, feature] = data[, paste0("AC_non_cancer", 
                                        ethnicityString, "_", gender[i])] - 
          data[, paste0("AC_controls_and_biobanks", 
                        ethnicityString, "_", gender[i])]
        
        feature = paste0("AN_non_cancer_non_control", ethnicityString, "_", 
                         gender[i])
        data[, feature] = data[, paste0("AN_non_cancer", 
                                        ethnicityString, "_", gender[i])] - 
          data[, paste0("AN_controls_and_biobanks", 
                        ethnicityString, "_", gender[i])]
      }
    }
    
    independentFeatures = grep(
  "^AC_cancer_.*X|^AC_non_cancer_non_control_.*X|^AC_controls_and_biobanks_.*X", 
      colnames(data), value = T)
    independentFeatures = sort(independentFeatures)
    independentFeaturesAN = sub("AC_", "AN_", independentFeatures)
  }
  
  # make sure all ANs are in the data
  if (sum(!(independentFeaturesAN %in% colnames(data))) > 0) {
    stop("AN IDs matched to AC IDs are not all in data")
  }
  
  dataIndep = t(data[, independentFeatures])
  ANIndep = round(rowMeans(t(data[, independentFeaturesAN])))
  
  # remove any variants with negative counts and print a warning
  countNegative = apply(dataIndep, MARGIN = 2, 
                        FUN = function(x) {any(is.na(x) | x < 0)})
  if (sum(countNegative) > 0) {
    indexNegative = which(countNegative)
    print("The following variants have negative counts (to be removed)")
    print(colnames(dataIndep)[indexNegative])
    if (length(indexNegative) == dim(dataIndep)[2]) {
      dataIndep = NULL
    } else {
      dataIndep = dataIndep[, -indexNegative]
    }
  }
  
  return(list(dataIndep = dataIndep, ANIndep = ANIndep))
}

identifyHighLDVariants = function(controlDataPerGene, 
                                  pThreshold = 0.05,
                                  ORThreshold = 1e3, 
                                  ACThreshold = 0,
                  ethnicitySet = c("nfe", "afr", "amr", "eas", "sas", "fin"),
                  gnomADVersion = "v2exome",
                  AFRatioThreshold = 1.5) {
  ## This function assumes using gnomAD as public controls
  
  # # <debug>
  # gene = "SLC26A10" #"RPL28" #"ZUFSP" #"CABYR" #"CFTR" #"CHIT1" #"MEOX2" #"KRTAP10-5" #"CEP192" #"KRT75" #"SLC26A10"
  # controlDataPerGene = controlCountQC[
  #    controlCountQC[, groupColumn] == gene & controlFunction, , drop = F]
  # rThreshold = 0.9
  # ratioThreshold = 0.8
  # pThreshold = 0.05
  # ORThreshold = 1e3
  # # </debug>
  
  data = NULL
  highLDList = NULL
  if (gnomADVersion == "v2genome") {
    ethnicitySet = setdiff(ethnicitySet, "sas")
  }
  # for (j in 1:length(ethnicitySet)) {
  #   ethnicity = ethnicitySet[j] 
  #   if (ethnicity != "") {
  #     ACString = paste0("AC_", ethnicity)
  #     ANString = paste0("AN_", ethnicity)
  #   } else {
  #     ACString = "AC"
  #     ANString = "AN"
  #   }
  #   features = grep(paste0("^", ACString, "|_", ACString, "|",
  #                          "^", ANString, "|_", ANString), 
  #                      colnames(controlDataPerGene), value = "T")
  #   if (is.null(data)) {
  #     data = controlDataPerGene[, features, drop = F]
  #     rownames(data) = controlDataPerGene[, 1]
  #   } else {
  #     data = cbind(data, controlDataPerGene[, features, drop = F])
  #   }
  #   
  #   # use uncorrelated counts
  #   gender = c("male", "female")
  #   for (i in 1:length(gender)) {
  #     # cancer 
  #     feature = paste0("cancer_", ACString, "_", gender[i])
  #     data[, feature] = data[, paste0(ACString, "_", gender[i])] - 
  #       data[, paste0("non_", feature)]
  #     
  #     feature = paste0("cancer_", ANString, "_", gender[i])
  #     data[, feature] = data[, paste0(ANString, "_", gender[i])] - 
  #       data[, paste0("non_", feature)]
  #     
  #     # non_control
  #     feature = paste0("non_cancer_non_control_", ACString, "_", gender[i])
  #     data[, feature] = 
  #       data[, paste0("non_cancer_", ACString, "_", gender[i])] - 
  #       data[, paste0("controls_", ACString, "_", gender[i])]
  #     
  #     feature = paste0("non_cancer_non_control_", ANString, "_", gender[i])
  #     data[, feature] = 
  #       data[, paste0("non_cancer_", ANString, "_", gender[i])] - 
  #       data[, paste0("controls_", ANString, "_", gender[i])]
  #   }
  # }
  #   
  # # independent counts
  # independentFeatures = grep(
  #   "^cancer_AC_.*male|non_cancer_non_control_AC_.*male|controls_AC_.*male", 
  #   colnames(data), value = T)
  # independentFeaturesAN = sub("_AC_", "_AN_", independentFeatures)
  # # make sure all ANs are in the data
  # if (sum(!(independentFeaturesAN %in% colnames(data))) > 0) {
  #   stop("AN IDs matched to AC IDs are not all in data")
  # }
  # 
  # dataIndep = t(data[, independentFeatures])
  # ANIndep = round(rowMeans(t(data[, independentFeaturesAN])))
  
  independentCounts = extractIndependentGnomADCounts(controlDataPerGene, 
                                                            ethnicitySet,
                                                     gnomADVersion)
  dataIndep = independentCounts[[1]]
  ANIndep = independentCounts[[2]]
  
  # filter and order by AC
  index2 = which(colSums(dataIndep) >= 2)
  if (length(index2) <= 1) {
    return(list(highLDList))
  }
  totalCounts = colSums(dataIndep[, index2], na.rm = T)
  # sort by the total counts of alleles, so the first variant can better 
  # represent the rest in high LD
  index2 = index2[order(totalCounts, decreasing = T)]
  dataUsed = dataIndep[, index2]
  ANUsed = ANIndep
  
  nVariants = dim(dataUsed)[2]
  satisfied = matrix(F, nVariants, nVariants)
  colnames(satisfied) = colnames(dataUsed)
  rownames(satisfied) = colnames(dataUsed)
  
  ##### use LRT based high LD detection ##########
  # calculate the log10(AF ratio)
  AFUsed = colSums(dataUsed) / sum(ANUsed)
  AFRatio = log10(as.matrix(AFUsed) %*% (1 / t(AFUsed)))
  AFRatio[lower.tri(AFRatio, diag = T)] = NA
  filteredList = which(abs(AFRatio) < log10(AFRatioThreshold), arr.ind = T)
  
  
  if (dim(filteredList)[1] == 0) {
    return(list(highLDList))
  }
  
  # estimate the OR between a pair of variant
  # # run on all pairs of variants
  # testResult = matrix(NA, nVariants * (nVariants - 1) / 2, 4)
  # index = 1
  # for (l in 1:(nVariants - 1)) {
  #   for (k in (l+1):nVariants) {
  #     dataForEstimation = cbind(dataUsed[, c(l, k)], n)
  #     result = LDTest(dataForEstimation, 3)
  #     testResult[index, ] = c(l, k, result[1:2])
  #     index = index + 1
  #   }
  # }

  variantID = colnames(dataUsed)
  # run on variants with similar AF
  if (dim(filteredList)[1] > 0) {
    for (i in 1:dim(filteredList)[1]) {
      l = filteredList[i, 1]
      k = filteredList[i, 2]
      dataForEstimation = cbind(dataUsed[, c(l, k)], ANUsed)
      result = LDTest(dataForEstimation, 
                                                  ORThreshold,
                                                  ACThreshold = ACThreshold)
      # # <debug>
      # if (result[1] > 0.1) {
      #   print(c(variantID[c(l, k)], result))
      # }
      # # </debug>
      if (result[1] < pThreshold) {
        satisfied[l, k] = T
        satisfied[k, l] = T
      }
    }
  }
    
  # extract connected components
  g  = graph.adjacency(satisfied)
  # plot(g)
  connectedComponents = components(g)
  if (sum(connectedComponents$csize > 1)) {
    cluster = which(connectedComponents$csize > 1)
    for (i in 1:length(cluster)) {
      highLDList[[i]] = names(which(connectedComponents$membership == 
                                      cluster[i]))
    }
  }
  
  return(list(highLDList))
}

XYNonPAR = function(chrID, position, reference) {
  # check whethere variants are within the non pseudoautosomal (non-PAR) region
  # of the X, Y chromosome

  # input
  # chrID: the chromosome ID of the variants
  # position: the position of the variants

  # output
  # a logical vector indicating whether each variant is in the non-PAR region

  if (reference == "GRCh37") {
    XPAR1 = c(60001, 2699520)
    YPAR1 = c(10001, 2699520)
    XPAR2 = c(154931044, 155260560)
    YPAR2 = c(59034050, 59363566)
  } else if (reference == "GRCh38") {
    XPAR1 = c(10001, 2781479)
    YPAR1 = c(10001, 2781479)
    XPAR2 = c(155701383, 156030895)
    YPAR2 = c(56887903,  57217415)
  } else {
    stop(paste("Reference build not supported:", reference))
  }

  inXYNonPAR = (grepl("X", chrID) & (
        position < XPAR1[1] | position > XPAR1[2] & position < XPAR2[1] |
                    position > XPAR2[2])) |
          (grepl("Y", chrID) & (
        position < YPAR1[1] | position > YPAR1[2] & position < YPAR2[1] | position > YPAR2[2])) 
  return(inXYNonPAR)
}