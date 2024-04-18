#' @name CoCoRV
#' @useDynLib CoCoRV, .registration = TRUE

#' @rawNamespace import(data.table, except = shift)
#' @importFrom S4Vectors Rle
#' @importFrom GenomicRanges GRanges strand reduce findOverlaps
#' @importFrom igraph components graph.adjacency make_undirected_graph
#' @import nloptr
#' @importFrom SeqArray seqOpen seqGetData seqApply seqSetFilter seqClose 
#'             seqSummary seqListVarData
#' @import Rcpp
#' @import RcppArmadillo
#' @import R.utils
#' @import IRanges
#' @import BiocParallel
#' @import DiscreteFDR
#' @import discreteMTP
NULL

#' Consistent summary Counts based Rare Variant burden test (CoCoRV)
#'
#' This is the main function for performing variant filtering, count extraction
#' or estimation, and the rare varaint burden test

#' @param controlGDSFile gds format file of controls
#' @param caseGDSFile gds format file of cases
#' @param controlAnnoGDSFile gds format file of annotated variants in controls
#' @param caseAnnoGDSFile gds format file of annotated variants in cases
#' @param nControl number of controls, will be detected based on AN if not set
#' @param nCase number of cases, will be detected based on AN if not set
#' @param sampleList a one column file listing sample IDs used in cases
#' @param outputPrefix the output file prefix
#' @param AFMax the maximum of the alternate alleles frequency
#' @param AFUse how to apply the AF threshold, options are: 
#'  case_control_joint (default), control, case_control_separate
#' @param AFGroup how to combine AFs from multiple groups, options are: 
#'  max: use the max of multiple group AFs (default), pool: pool all groups 
#'  together to calculate the AF
#' @param countUse how to use the count for association test. 
#'  sample: each sample is a unit in the final association test (default), 
#'  or allele_null: functional allele counts vs nonfunctional allele counts as 
#'  in proxecat
#' @param variantMissing the maximum missingness allowed for a variant
#' @param removeStar remove overlapping deletion variants with the star symbol.
#'  This is useful when using gnomAD because gnomAD does not have variants with
#'  star as alternate allele.
#' @param groupColumn the column specifying the column name grouping variants, 
#' such as the annotated gene name. If NA, all variants will be considered.
#' @param ACANConfig a file specifying the AC AN info ID and optional 
#' stratified ID for cases
#' @param caseGroup a file with sample ID in the first column and the
#' stratified group info in the 2nd column
#' @param variantExcludeFile a one column file specifying the variants to be 
#' excluded
#' @param variantIncludeFile a one column file specifying the variants to be 
#' included, it will also include variants specified in variantExcludeFile
#' @param ignoreFilter if not specified, remove variants with filter not equal 
#' to PASS in either case or control for a consistent filtering, default is 
#' false.  
#' @param variantGroupCustom an optional R file with functions defining the 
#' variants of interest or a tab separated two column text file where the 
#' first column is the variantGroup ID and the second column is the filtering
#' expression using annotation column names. 
#' @param extraParamJason an optional JASON file to provide extra parameters
#' to the custuomized file variantGroupCustom
#' @param variantGroup a specified variant group to use for test or a self 
#' defined function to define the variants of interest. See below for more 
#' details on defined variant sets. 
#' @param annotationUsed a one column file without the header listing 
#' all required annotations used in variant filtering. AC AN annotations from 
#' ACANConfig will be added automatically.
#' @param minREVEL the minimum REVEL score for pathogenic missense variants
#' @param maxAFPopmax the maximum AF_popmax annotated from gnomAD, make sure to 
#' use the same annotated AF_popmax for both the case and the control. 
#' To do this, first remove the annotation AF_popmax from the raw gnomAD vcf 
#' using bcftools annotate and then add the annotation from the 
#' annotation tool, e.g., ANNOVAR 
#' @param bed a bed file to specify the regions considered, e.g., the 
#' intersection of the bed files from cases and controls, where >= 90% of 
#' samples have coverage >= 10
#' @param overlapType overlap type between the variant and the bed file: 
#' "any", "within" (default), "gds"
#' @param fullCaseGenotype If true (default), use the full genotypes to extract 
#' counts for cases, otherwise summary counts are used
#' @param checkHighLDInControl for each case double het, if high LD in control 
#' then exclude
#' @param gnomADVersion specify the major gnomAD version, v2exome (default),
#' v2genome, v3genome,
#'  if gnomAD is used
#' @param pLDControl the pvalue threshold to detect high LD variants in control
#' @param ORThresholdControl the OR threshold to define high LD variants
#' @param highLDVariantFile the precomputed variant pairs in LD in controls
#' @param ignoreEthnicityInLD whether to ignore ethnicity when loading the 
#' precomputed LD pairs 
#' @param regionFile a bed file with names at the 4th column defining the 
#' processing unit. This can be used 
#' to define test regions or to process relatively small regions 
#' sequentially to reduce memory
#' @param reference the reference build, either GRCh37 or GRCh38, used for 
#'  defining the pseudoautosomal region in X and Y

#' @return A list of three components
#'           testResult: the association test result, 
#'           caseGeneList: the gene and variants used in the test in cases
#'           controlGeneList: the gene and variants used in the test in controls

#' @details 
#' Defined variant sets based on annovar annotations
#' \itemize{
#'   \item annovar_pathogenic: frameshift_deletion | frameshift_insertion | 
#'     stopgain | nonsynonymous_SNV with REVEL >= minREVEL   
#'   \item annovar_function: frameshift_deletion | frameshift_insertion | 
#'     stopgain | splicing | nonsynonymous_SNV with REVEL >= minREVEL
#'   \item annovar_missense: nonsynonymous_SNV with REVEL >= minREVEL
#'   \item annovar_LOF: frameshift_deletion | frameshift_insertion | stopgain
#'   \item annovar_synonym: synonymous_SNV
#'   \item annovar_splicing: splicing
#'   \item annovar_intron: intronic | ncRNA_intronic
#'   \item annovar_UTR3: UTR3
#'   \item annovar_UTR5: UTR5
#'   \item annovar_UTR: UTR3 | UTR5  
#' }
#' For a self defined function, it must take a data frame of summary count data
#' as the parameter and return a logical index on rows, see customVariantSets
#' in variantSets.R for an example. 
#'
#' If stratified analysis is specified, the AC and AN will still be used for
#' filtering and QC, then the counts corresponding to the stratified groups
#' will be used in the CMH exact test

#' @export
CoCoRV = function(controlGDSFile,
                  caseGDSFile,
                  controlAnnoGDSFile = "",
                  caseAnnoGDSFile = "",
                  nControl = -1,
                  nCase = -1,
                  sampleList = "",
                  outputPrefix = "out",
                  AFMax = 1e-3,
                  AFUse = "case_control_joint",
                  AFGroup = "max",
                  countUse = "sample",
                  variantMissing = 0.1,
                  removeStar = F,
                  minREVEL = 0.65,
                  maxAFPopmax = 1,
                  ACANConfig = "",
                  caseGroup = "",
                  variantExcludeFile = "",
                  variantIncludeFile = "",
                  ignoreFilter = F,
                  variantGroupCustom = "",
                  extraParamJason = "",
                  variantGroup = "pathogenic",
                  groupColumn = "Gene.refGene",
                  annotationUsed = "",
                  bed = "",
                  fullCaseGenotype = T,
                  overlapType = "within",
                  checkHighLDInControl = F,
                  gnomADVersion = "v2exome",
                  pLDControl = 0.05,
                  ORThresholdControl = 1,
                  highLDVariantFile = "",
                  ignoreEthnicityInLD = T,
                  regionFile = "",
                  reference = "GRCh37") {
  variantExclude = NULL
  if (nchar(variantExcludeFile) > 0) {
    # variantExclude = c(as.matrix(read.table(variantExcludeFile)))
    variantExclude = fread(variantExcludeFile,data.table = F, header = F, 
                           colClasses="character")[, 1]
  }

  variantInclude = NULL
  if (nchar(variantIncludeFile) > 0) {
    variantInclude = fread(variantIncludeFile,data.table = F, header = F, 
                           colClasses="character")[, 1]
  }
  
  sampleID = NULL
  if (sampleList != "") {
    sampleID = c(as.matrix(read.table(sampleList, header = F, as.is = T)))
  }
  
  bedGRange = NULL
  if (bed != "") {
    #bedRegions = read.table(bed, header = F, as.is = T)
    bedRegions = fread(bed, na.strings = ".", data.table = F)
    bedGRange = GRanges(seqnames = bedRegions[, 1],
                        ranges = IRanges(start = bedRegions[, 2] + 1,
                                         end = bedRegions[, 3]),
                        strand = Rle(strand("+"), dim(bedRegions)[1]))
    bedGRange = reduce(bedGRange)
  }
  
  regionGRange = NULL
  if (regionFile != "") {
    region = fread(regionFile, na.strings = ".", 
                   data.table = F)
    if (dim(region)[1] >= 1) {
      regionGRange = GRanges(seqnames = region[, 1],
                             ranges = IRanges(start = region[, 2] + 1,
                                              end = region[, 3]),
                             strand = Rle(strand("+"), dim(region)[1]))
    }
  } 
  
  ACANConfigData = NULL
  if (ACANConfig != "") {
    ACANConfigData = as.matrix(read.table(ACANConfig, header = T, as.is = T))
    if ("caseGroup" %in% colnames(ACANConfigData)) {
      # make sure the caseGroup is unqiue 
      stopifnot(dim(ACANConfigData)[1] == 
                  length(unique(ACANConfigData[, "caseGroup"])))
    }
  }
  
  caseGroupInfo = NULL
  sexID = NULL
  if (caseGroup != "" && ("caseGroup" %in% colnames(ACANConfigData))) {
    caseInfo = as.matrix(read.table(caseGroup, header = T, as.is = T))
    indexUsed = match(sampleID, caseInfo[, 1])
    if (sum(is.na(indexUsed)) > 0) {
      stop("some sample IDs do not have matched group info")
    }
    caseGroupInfo = caseInfo[indexUsed, 2]
    if ("SEX" %in% colnames(caseInfo)) {
      # generate sexID matched to sampleID
      sexID = caseInfo[indexUsed, "SEX"]
    }
  }
  
  if (caseGroup == "" && fullCaseGenotype) {
    if (!is.null(sampleID)) {
      # assume all case samples are from the same group
      stopifnot(dim(ACANConfigData)[1] == 1)
      caseGroupInfo = rep("default_group", length(sampleID))
      ACANConfigData = cbind(ACANConfigData, "default_group")
      colnames(ACANConfigData)[dim(ACANConfigData)[2]] = "caseGroup"
    } else {
      # need to make sure caseAC caseAN are specified in the config file
      stopifnot(c("caseAC", "caseAN") %in% colnames(ACANConfigData))
    }
  }

  # check annotationUsed
  # predefined groups
  definedSet = c("annovar_pathogenic", 
                 "annovar_function", 
                 "annovar_LOF", 
                 "annovar_missense",
                 "annovar_synonym",
                 "annovar_splicing", 
                 "annovar_intron", 
                 "annovar_UTR3",
                 "annovar_UTR5", 
                 "annovar_UTR")
  if (annotationUsed == "") {
    if (variantGroup %in% definedSet) {
      annotationSet = c("ExonicFunc.refGene", "Func.refGene", "REVEL")
      if (groupColumn != "NA") {
        annotationSet = union(annotationSet, groupColumn)
      }
    } else {
      cat("!! All annotations from data will be loaded. It is much faster and more memory efficient to specify needed annotations by using the option annotationUsed\n")
      annotationSet = NULL
    }
  } else {
    annotationSet = c(as.matrix(
             read.table(annotationUsed, header = F, as.is = T)))
    if (groupColumn != "NA") {
        annotationSet = union(annotationSet, groupColumn)
    }
  }
  
  
  # open gds files
  controlGDS = NULL
  caseGDS = NULL
  if (controlGDSFile != "") {
    controlGDS = seqOpen(controlGDSFile)
  }
  if (caseGDSFile != "") {
    caseGDS = seqOpen(caseGDSFile)
  }

  controlAnnoGDS = NULL
  caseAnnoGDS = NULL
  if (controlAnnoGDSFile != "") {
    if (controlAnnoGDSFile == controlGDSFile) {
      controlAnnoGDS = controlGDS
    } else {
      controlAnnoGDS = seqOpen(controlAnnoGDSFile)
    }
  }
  if (caseAnnoGDSFile != "") {
    if (caseAnnoGDSFile == caseGDSFile) {
      caseAnnoGDS = caseGDS
    } else {
      caseAnnoGDS = seqOpen(caseAnnoGDSFile)
    }
  }
  
  highLDLoaded = NULL
  if (highLDVariantFile != "") {
    cat("load precomputed high LD variants\n")
    highLDLoaded = loadHighLDVariants(highLDVariantFile, ignoreEthnicityInLD, 
                                      pThreshold = pLDControl,
                                      ORThreshold = ORThresholdControl)
    cat("finished loading of high LD variants\n")
  }
  
  resultPerRegion = list()
  k = 1
  if (!is.null(regionGRange)) {
    for (i in 1:length(regionGRange)) {
      bedPerRegion = intersect(regionGRange[i], bedGRange)
      if (length(bedPerRegion) > 0) {
        result = CoCoRVPerRegion(controlGDS,
                              caseGDS,
                              controlAnnoGDS = controlAnnoGDS,
                              caseAnnoGDS = caseAnnoGDS,
                              nControl = nControl,
                              nCase = nCase,
                              sampleID = sampleID,
                              outputPrefix = outputPrefix,
                              AFMax = AFMax,
                              AFUse = AFUse,
                              AFGroup = AFGroup,
                              countUse = countUse,
                              variantMissing = variantMissing,
                              removeStar = removeStar,
                              minREVEL = minREVEL,
                              maxAFPopmax = maxAFPopmax,
                              ACANConfigData = ACANConfigData,
                              caseGroupInfo = caseGroupInfo,
                              variantExclude = variantExclude,
                              variantInclude = variantInclude,
                              ignoreFilter = ignoreFilter,
                              variantGroupCustom = variantGroupCustom,
                              variantGroup = variantGroup,
                              extraParamJason=extraParamJason,
                              groupColumn = groupColumn,
                              annotationSet = annotationSet,
                              bedGRange = bedPerRegion,
                              fullCaseGenotype = fullCaseGenotype,
                              overlapType = overlapType,
                              checkHighLDInControl = checkHighLDInControl,
                              gnomADVersion = gnomADVersion,
                              pLDControl = pLDControl,
                              ORThresholdControl = ORThresholdControl,
                              highLDVariantFile = highLDVariantFile,
                              ignoreEthnicityInLD = ignoreEthnicityInLD,
                              regionGRange = regionGRange[i],
                              highLDLoaded = highLDLoaded,
                              reference = reference,
                              sexID = sexID)
        result$regionID = region[i, 4]
        resultPerRegion[[k]] = result
        k = k + 1
        gc()
      }
    }
  } else {
    resultPerRegion[[1]] = CoCoRVPerRegion(controlGDS,
                             caseGDS,
                             controlAnnoGDS = controlAnnoGDS,
                             caseAnnoGDS = caseAnnoGDS,
                             nControl = nControl,
                             nCase = nCase,
                             sampleID = sampleID,
                             outputPrefix = outputPrefix,
                             AFMax = AFMax,
                             AFUse = AFUse,
                             AFGroup = AFGroup,
                             countUse = countUse,
                             variantMissing = variantMissing,
                             removeStar = removeStar,
                             minREVEL = minREVEL,
                             maxAFPopmax = maxAFPopmax,
                             ACANConfigData = ACANConfigData,
                             caseGroupInfo = caseGroupInfo,
                             variantExclude = variantExclude,
                             variantInclude = variantInclude,
                             ignoreFilter = ignoreFilter,
                             variantGroupCustom = variantGroupCustom,
                             variantGroup = variantGroup,
                             extraParamJason = extraParamJason,
                             groupColumn = groupColumn,
                             annotationSet = annotationSet,
                             bedGRange = bedGRange,
                             fullCaseGenotype = fullCaseGenotype,
                             overlapType = overlapType,
                             checkHighLDInControl = checkHighLDInControl,
                             gnomADVersion = gnomADVersion,
                             pLDControl = pLDControl,
                             ORThresholdControl = ORThresholdControl,
                             highLDVariantFile = highLDVariantFile,
                             ignoreEthnicityInLD = ignoreEthnicityInLD,
                             regionGRange = NULL,
                             highLDLoaded = highLDLoaded,
                             reference = reference,
                             sexID = sexID)
  }
  
  # close gds file
  if (!is.null(caseGDS)) {
    seqClose(caseGDS)
  }
  if (!is.null(controlGDS)) {
    seqClose(controlGDS)
  }
  if (!is.null(caseAnnoGDS) && caseAnnoGDSFile != caseGDSFile) {
    seqClose(caseAnnoGDS)
  }
  if (!is.null(controlAnnoGDS) && controlAnnoGDSFile != controlGDSFile) {
    seqClose(controlAnnoGDS)
  }

  return(resultPerRegion)
}

CoCoRVPerRegion = function(controlGDS,
                  caseGDS,
                  controlAnnoGDS = NULL,
                  caseAnnoGDS = NULL,
                  nControl = -1,
                  nCase = -1,
                  sampleID = NULL,
                  outputPrefix = "out",
                  AFMax = 1e-3,
                  AFUse = "case_control_joint",
                  AFGroup = "max",
                  countUse = "sample",
                  variantMissing = 0.1,
                  removeStar = F,
                  minREVEL = 0.65,
                  maxAFPopmax = 1,
                  ACANConfigData = NULL,
                  caseGroupInfo = NULL,
                  variantExclude = NULL,
                  variantInclude = NULL,
                  ignoreFilter = F,
                  variantGroupCustom = "",
                  variantGroup = "pathogenic",
                  extraParamJason = "",
                  groupColumn = "Gene.refGene",
                  annotationSet = NULL,
                  bedGRange = NULL,
                  fullCaseGenotype = T,
                  overlapType = "within",
                  checkHighLDInControl = F,
                  gnomADVersion = "v2exome",
                  pLDControl = 0.05,
                  ORThresholdControl = 1,
                  highLDVariantFile = "",
                  ignoreEthnicityInLD = F,
                  regionGRange = NULL,
                  highLDLoaded = NULL,
                  reference = "GRCh37",
                  sexID = NULL) {
  
  # There are two ways to specify the configurations of the analysis
  # 1. ACANConfig & caseGroupInfo
  #   ACANConfig includes column headers controlAC controlAN caseGroup which 
  #   are the AC AN Info ID in controls and corresponding case group ID. 
  #   caseGroupInfo stores the case group IDs for each sample. 
  #   This configuration requires full genotypes of the cases
  # 2. ACANConfig only
  #   This includes column headers controlAC controlAN caseAC caseAN
  #   This does not require the full genotypes of cases, support only one line
  #   in the configuration file for now. This does not require full genotypes
  #   of the cases
  stopifnot((!is.null(ACANConfigData) && !is.null(caseGroupInfo)) || 
             (!is.null(ACANConfigData) && 
                ("caseAC" %in% colnames(ACANConfigData) && 
                ("caseAN" %in% colnames(ACANConfigData)))
             ) || 
             (is.null(ACANConfigData) && is.null(caseGroupInfo))
           )
  
  caseAC = "AC"
  caseAN = "AN"
  controlAC = "AC"
  controlAN = "AN"
  caseAA = NULL
  controlAA = NULL
  if (!is.null(ACANConfigData) && ("caseAC" %in% colnames(ACANConfigData)) && 
      ("caseAN" %in% colnames(ACANConfigData))) {
    caseAC = ACANConfigData[, "caseAC"]
    caseAN = ACANConfigData[, "caseAN"]
    caseAA = NULL
    if ("caseAA" %in% colnames(ACANConfigData)) {
      caseAA = ACANConfigData[, "caseAA"]
    }
    controlAC = ACANConfigData[, "controlAC"]
    controlAN = ACANConfigData[, "controlAN"]
    controlAA = NULL
    if ("controlAA" %in% colnames(ACANConfigData)) {
      controlAA = ACANConfigData[, "controlAA"]
    }
  } else if (!is.null(ACANConfigData)) {
    # only controls AC AN are specified
    controlAC = ACANConfigData[, "controlAC"]
    controlAN = ACANConfigData[, "controlAN"]
    controlAA = NULL
    if ("controlAA" %in% colnames(ACANConfigData)) {
      controlAA = ACANConfigData[, "controlAA"]
    }
    # match the number of groups in cases
    caseAC = rep(caseAC, length(controlAC))
    caseAN = rep(caseAN, length(controlAC))
  }
  
  caseCountFile = ""
  controlCountFile = ""
  if (AFUse == "control") {
    AFControlOnly = T 
  } else if (AFUse == "case_control_joint" || 
             AFUse == "case_control_separate") {
    AFControlOnly = F
  } else {
    stop("wrong AFUse values. 
         Can only be case_control_joint, control, case_control_separate")
  }
  
  # extract data
  if (!is.null(controlGDS)) {
    ACANIDs = c(controlAC, controlAN, controlAA)
    if (is.null(controlAnnoGDS)) {
      annotationSetFromControlGDS = annotationSet
    } else {
      annotationSetFromControlGDS = ACANIDs
    }
    cat("loading control count data ...\n")
    controlCount = extractACANAnnotations(controlGDS, sampleID = NULL, 
                                      bedGRange, 
                                      updateACAN = F,
                                      overlapType = overlapType, 
                                      ACANIDs = ACANIDs, 
                                  annotationSet = annotationSetFromControlGDS)
  } else {
    controlCount = fread(controlCountFile, sep = "\t", na.strings = ".",
                         data.table = F)
  }
  if (is.null(controlCount)) {
    return(NULL)
  }


  if (!is.null(caseGDS) && caseCountFile == "") {
    if (fullCaseGenotype) {
      if (is.null(caseAnnoGDS)) {
        annotationSetFromCaseGDS = annotationSet
      } else {
        annotationSetFromCaseGDS = c("AC", "AN")
      }
      cat("loading case count data ...\n")
      caseCount = extractACANAnnotations(caseGDS, sampleID, 
                                         bedGRange,
                                       updateACAN = T, 
                                       overlapType = overlapType,
                                    annotationSet = annotationSetFromCaseGDS)
    } else {
      ACANIDs = c(caseAC, caseAN, caseAA)
      if (is.null(caseAnnoGDS)) {
        annotationSetFromCaseGDS = annotationSet
      } else {
        annotationSetFromCaseGDS = ACANIDs
      }
      caseCount = extractACANAnnotations(caseGDS, sampleID = NULL, 
                                         bedGRange,
                                         updateACAN = F, 
                                         overlapType = overlapType,
                                         ACANIDs = ACANIDs, 
                                   annotationSet = annotationSetFromCaseGDS)
    }
  } else {
    caseCount = fread(caseCountFile, sep = "\t", na.strings = ".", 
                      data.table = F) 
  }
  if (is.null(caseCount)) {
    return(NULL)
  }

  # for pure variant annotations only, no AC AN counts loaded 
  caseAnno = NULL
  if (!is.null(caseAnnoGDS)) {
    cat("loading case annotation data ...\n")
    caseAnno = extractACANAnnotations(caseAnnoGDS, sampleID, 
                                       bedGRange,
                                       updateACAN = F, 
                                       overlapType = overlapType,
                                       ACANIDs = NULL, 
                                       annotationSet = annotationSet) 
  }
  controlAnno = NULL
  if (!is.null(controlAnnoGDS)) {
    cat("loading control annotation data ...\n")
    controlAnno = extractACANAnnotations(controlAnnoGDS, sampleID = NULL, 
                                      bedGRange, 
                                      updateACAN = F,
                                      overlapType = overlapType,
                                      ACANIDs = NULL, 
                                      annotationSet = annotationSet)
  }

  extractStratifiedCounts = !is.null(ACANConfigData) & 
    !is.null(caseGroupInfo) & 
    !(("caseAC" %in% colnames(ACANConfigData)) & 
    ("caseAN" %in% colnames(ACANConfigData))) & 
    (length(controlAC) >= 1) & 
    fullCaseGenotype
  
  caseACANExtracted = NULL
  nVariant = dim(caseCount)[1]
  if (extractStratifiedCounts) {
    groupIDSpecified = ACANConfigData[, "caseGroup"]
    caseACANExtracted = extractACANCounts(caseGDS, sampleID, 
                                          bedGRange,
                                caseGroupInfo, groupIDSpecified, nVariant,
                                overlapType = overlapType, 
                                sexID = sexID,
                                reference = reference)
  }
  
  # set NA to 0 for AC AN counts
  # controls
  index = which(grepl("AC$|AN$|_AC_|_AN_|^AC_|^AN_|nhomalt", 
                      colnames(controlCount)) & 
                  !grepl("popmax", colnames(controlCount)))
  if (sum(is.na(controlCount[, index])) > 0) {
    print("NA exist in control counts, set NA to 0 but better confirm")
    for (i in 1:length(index)) {
      indexNA = is.na(controlCount[, index[i]])
      if (sum(indexNA) > 0) {
        controlCount[indexNA, index[i]] = 0
      }
    }
  }
  # cases
  index = which(grepl("AC$|AN$|_AC_|_AN_|^AC_|^AN_|nhomalt", 
                      colnames(caseCount)) & 
                  !grepl("popmax", colnames(caseCount)))
  if (sum(is.na(caseCount[, index])) > 0) {
    print("NA exist in control counts, set NA to 0 but better confirm")
    for (i in 1:length(index)) {
      indexNA = is.na(caseCount[, index[i]])
      if (sum(indexNA) > 0) {
        caseCount[indexNA, index[i]] = 0
      }
    }
  }
  
  for (i in 1:length(controlAC)) {
    controlCount[, controlAC[i]] = as.numeric(controlCount[, controlAC[i]])
    controlCount[, controlAN[i]] = as.numeric(controlCount[, controlAN[i]])
    caseCount[, caseAC[i]] = as.numeric(caseCount[, caseAC[i]])
    caseCount[, caseAN[i]] = as.numeric(caseCount[, caseAN[i]])
  }
  
  # convert AF_popmax to double
  if ("AF_popmax" %in% colnames(controlCount)) {
    controlCount[, "AF_popmax"] = as.numeric(controlCount[, "AF_popmax"])
  }
  if ("AF_popmax" %in% colnames(caseCount)) {
    caseCount[, "AF_popmax"] = as.numeric(caseCount[, "AF_popmax"])
  }
  
  caseCountQC = caseCount
  controlCountQC = controlCount
  
  # estimate AF used for filtering
  if (length(controlAC) == 1) {  # one group
    if (AFUse == "case_control_joint") {
      # estimate AF based on both cases and controls instead 
      estimatedAF = estimateAFJoint(controlCountQC[, 1], 
                                    controlCountQC[, controlAC],
                                    controlCountQC[, controlAN],
                                    caseCountQC[, 1], caseCountQC[, caseAC], 
                                    caseCountQC[, caseAN])
      controlAF = estimatedAF$controlAF
      caseAF = estimatedAF$caseAF
    } else {
      # separate estimation of AF
      controlAF = controlCountQC[, controlAC] / controlCountQC[, controlAN]
      caseAF = caseCountQC[, caseAC] / caseCountQC[, caseAN]
    }
  } else {
    # multiple groups, only support case_control_joint
    if(AFUse != "case_control_joint") {
      stop("only joint AF estimate is supported for multiple groups")
    }
    if (AFGroup == "pool") {
      controlCountQCPooledAC = integer(dim(controlCountQC)[1])
      controlCountQCPooledAN = integer(dim(controlCountQC)[1])
      caseCountQCPooledAC = integer(dim(caseCountQC)[1])
      caseCountQCPooledAN = integer(dim(caseCountQC)[1])
      for (i in 1:length(controlAC)) {
        controlCountQCPooledAC = controlCountQCPooledAC + 
                                   controlCountQC[, controlAC[i]]
        controlCountQCPooledAN = controlCountQCPooledAN + 
          controlCountQC[, controlAN[i]]
        if (extractStratifiedCounts) {
          caseCountQCPooledAC = caseCountQCPooledAC + 
            caseACANExtracted[[1]][, i]
          caseCountQCPooledAN = caseCountQCPooledAN + 
            caseACANExtracted[[2]][, i]
        } else {
          caseCountQCPooledAC = caseCountQCPooledAC + 
            caseCountQC[, caseAC[i]]
          caseCountQCPooledAN = caseCountQCPooledAN + 
            caseCountQC[, caseAN[i]]
        }
      }
      estimatedAF = estimateAFJoint(controlCountQC[, 1], 
                                    controlCountQCPooledAC,
                                    controlCountQCPooledAN,
                                    caseCountQC[, 1], caseCountQCPooledAC, 
                                    caseCountQCPooledAN)
      controlAF = estimatedAF$controlAF
      caseAF = estimatedAF$caseAF
    } else if (AFGroup == "max") {
      controlAFMatrix = matrix(0L, dim(controlCountQC)[1], length(controlAC))
      caseAFMatrix = matrix(0L, dim(caseCountQC)[1], length(controlAC))
      for (i in 1:length(controlAC)) {
        if (!extractStratifiedCounts) {
          caseACUsed = caseCountQC[, caseAC[i]]
          caseANUsed = caseCountQC[, caseAN[i]]
        } else {
          caseACUsed = caseACANExtracted[[1]][, i]
          caseANUsed = caseACANExtracted[[2]][, i]
        }
        estimatedAF = estimateAFJoint(controlCountQC[, 1], 
                                      controlCountQC[, controlAC[i]],
                                      controlCountQC[, controlAN[i]],
                                      caseCountQC[, 1], 
                                      caseACUsed, 
                                      caseANUsed)
        controlAFMatrix[, i] = estimatedAF$controlAF
        caseAFMatrix[, i] = estimatedAF$caseAF
      }
      
      controlAF = apply(controlAFMatrix, MARGIN = 1, max)
      caseAF = apply(caseAFMatrix, MARGIN = 1, max)
    } else {
      stop("AFGroup must be either max or pool")
    }
  }
  
  # filter by missingness in the controls and cases
  # pool all ANs together if there are multiple groups
  # This should be done first before any filtering because it relies on 
  # the matching of caseACANExtracted and caseCountQC
  controlANMatrix = controlCountQC[, controlAN, drop = F]
  controlANs = apply(controlANMatrix, MARGIN = 1, FUN = sum)
  # handle nonPAR regions in the sex chromosomes
  controlVariantID = controlCountQC[, 1]
  variantInfo = matrix(unlist(strsplit(controlVariantID, "-")), 
      ncol = 4, byrow = T)
  chrID = variantInfo[, 1]
  position = as.integer(variantInfo[, 2])
  variantInNonPAR = XYNonPAR(chrID, position, reference)
  controlANMax = numeric(length(controlANs))
  controlANMax[] = max(controlANs)
  if (sum(variantInNonPAR) > 0) {
    controlANMax[variantInNonPAR] = max(controlANs[variantInNonPAR])
  }
  

  controlRemove = (controlANs <  controlANMax * (1 - variantMissing))
  IDRemoveInControl = controlCountQC[controlRemove, 1]
  if (!extractStratifiedCounts) {
    caseANMatrix = caseCountQC[, caseAN, drop = F]
  } else {
    caseANMatrix = caseACANExtracted[[2]]
  }
  caseANs = apply(caseANMatrix, MARGIN = 1, FUN = sum)
  caseVariantID = caseCountQC[, 1]
  variantInfo = matrix(unlist(strsplit(caseVariantID, "-")), 
      ncol = 4, byrow = T)
  chrID = variantInfo[, 1]
  position = as.integer(variantInfo[, 2])
  variantInNonPAR = XYNonPAR(chrID, position, reference)
  caseANMax = numeric(length(caseANs))
  caseANMax[] = max(caseANs)
  if (sum(variantInNonPAR) > 0) {
    caseANMax[variantInNonPAR] = max(caseANs[variantInNonPAR])
  }

  caseRemove = (caseANs < caseANMax * (1 - variantMissing))
  IDRemoveInCase = caseCountQC[caseRemove, 1]
  IDRemove = union(IDRemoveInControl, IDRemoveInCase)
  
  variantExclude = union(variantExclude, IDRemove)
  
  # filter by AF
  ########## filter on controls ###########
  indexControlRemove = (controlAF > AFMax) & 
                 !(controlCountQC[, 1] %in% variantInclude)
  variantIDRemove = controlCountQC[indexControlRemove, 1]
  controlCountQC = controlCountQC[!indexControlRemove, , drop = F]
  ####### filter on cases ###########
  if (AFControlOnly) {
    indexCaseRemove = (caseCountQC[, 1] %in% variantIDRemove)
  } else {
    indexCaseRemove = (caseAF > AFMax | 
                         caseCountQC[, 1] %in% variantIDRemove)
  }
  indexCaseRemove = indexCaseRemove & !(caseCountQC[, 1] %in% variantInclude)
  caseCountQC = caseCountQC[!indexCaseRemove, , drop = F]
  
  # remove SNPs with * 
  if (removeStar) {
    # remove star in cases
    caseRemove = grep("-*", caseCountQC[, 1], fixed = T, value = T)
    # remove star in controls
    controlRemove = grep("-*", controlCountQC[, 1], fixed = T, value = T)
    IDRemove = union(caseRemove, controlRemove)
    variantExclude = union(variantExclude, IDRemove)
  }
  
  # This step excludes any variants that failed in any cohort
  if (!ignoreFilter) {  
    if (!("FILTER" %in% colnames(caseCountQC))) {
      stop("FILTER column is not in extracted annotations")
    }
    controlFailed = controlCountQC[controlCountQC$FILTER != "PASS", 1]
    caseFailed = caseCountQC[caseCountQC$FILTER != "PASS", 1]
    variantFailed = union(controlFailed, caseFailed)
    variantExclude = union(variantExclude, variantFailed)
  }

  variantExclude = setdiff(variantExclude, variantInclude)
  
  # remove variants in variantExclude
  if (length(variantExclude) > 0) {
    indexRemoveControl = (controlCountQC[, 1] %in% variantExclude)
    controlCountQC = controlCountQC[!indexRemoveControl, , drop = F]
    indexRemoveCase = (caseCountQC[, 1] %in% variantExclude)
    caseCountQC = caseCountQC[!indexRemoveCase, , drop = F]
  }
  
  # remove duplicated
  indexDuplicated = duplicated(controlCountQC[, 1])
  if (sum(indexDuplicated) > 0) {
    controlCountQC = controlCountQC[!duplicated(controlCountQC[, 1]), , 
                                    drop = F]
  }
  indexDuplicated = duplicated(caseCountQC[, 1])
  if (sum(indexDuplicated) > 0) {
    caseCountQC = caseCountQC[!indexDuplicated, , drop = F]
  }

  if (!is.null(caseAnno)) {
    index = match(caseCountQC[, 1], caseAnno[, 1])
    if (sum(is.na(index)) > 0) {
      stop("annotation data do not cover all variants in count data for cases")
    }
    caseAnnotation = caseAnno[index, , drop = F]

    # remove annotations existing in caseAnnotation from caseCountQC 
    annoDuplicated = (colnames(caseCountQC) %in% colnames(caseAnnotation))
    annoDuplicated[1] = F # keep the variant ID
    caseCountQC = caseCountQC[, !annoDuplicated]
    # add annotation to caseCountQC
    caseCountQC = cbind(caseCountQC, caseAnnotation[, -1, drop = F])
  } 
  if (!is.null(controlAnno)) {
    index = match(controlCountQC[, 1], controlAnno[, 1])
    if (sum(is.na(index)) > 0) {
      stop("annotation data do not cover all variants in count data for controls")
    }
    controlAnnotation = controlAnno[index, , drop = F]

    # remove annotations existing in controlAnnotation from controlCountQC 
    annoDuplicated = (colnames(controlCountQC) %in% 
                         colnames(controlAnnotation))
    annoDuplicated[1] = F # keep the variant ID
    controlCountQC = controlCountQC[, !annoDuplicated]
    # add annotation to caseCountQC
    controlCountQC = cbind(controlCountQC, controlAnnotation[, -1, drop = F])
  } 
  
  # define functional variants
  controlFunction = variantSets(controlCountQC, extraParamJason, 
                                variantGroupCustom,
                                variantGroup, minREVEL, maxAFPopmax) | 
                    (controlCountQC[, 1] %in% variantInclude)
  caseFunction = variantSets(caseCountQC, extraParamJason, 
                             variantGroupCustom,
                             variantGroup, minREVEL, maxAFPopmax) | 
                    (caseCountQC[, 1] %in% variantInclude)
  if (sum(controlFunction) == 0 && sum(caseFunction) == 0) {
    stop("no variants selected in either case or control")
  }
  
  complementary = F
  benignUpperBound = 0.2
  caseNonFunction = NULL
  controlNonFunction = NULL
  if (regexpr("^sample", countUse) == -1) {
    if (complementary) { 
      # set all the rest to nonFunction
      caseNonFunction = !caseFunction
      controlNonFunction = !controlFunction
    } else {
      functional = "ExonicFunc.refGene"
      nonfunctionSet = "synonymous_SNV"
      pathogenicScore = "REVEL"
      caseNonFunction = (caseCountQC[, functional] %in% nonfunctionSet | 
                           (!is.na(caseCountQC[, pathogenicScore]) & 
                        caseCountQC[, pathogenicScore] < benignUpperBound)
      )
      controlNonFunction = (controlCountQC[, functional] %in% nonfunctionSet | 
                              (!is.na(controlCountQC[, pathogenicScore]) & 
                         controlCountQC[, pathogenicScore] < benignUpperBound))
    }
  }

  ## process each gene
  if (groupColumn == "NA") {  
    print("groupColumn is NA, all variants will be considered")
    groupColumn = "__tempGroup__"
    controlCountQC[, "__tempGroup__"] = "tempGroup"
    caseCountQC[, "__tempGroup__"] = "tempGroup"
  }
  geneList = union(unique(controlCountQC[, groupColumn]), 
                   unique(caseCountQC[, groupColumn]))
  nGene = length(geneList)
  controlACANConfig = NULL
  caseACANConfig = NULL
  nControlPerGroup = NULL
  nCasePerGroup = NULL
  nGroup = NULL
  
  if (nCase == -1 || nControl == -1) {
    print("Number of cases and controls in each group")
    # assume there are always other regions besides XYNonPAR
    # TODO: may need to handle the situations where there are only XYNonPar
    #   regions, even though this is likely to be rare
    variantInfo = matrix(unlist(strsplit(controlCountQC[, 1], "-")), 
      ncol = 4, byrow = T)
    chrID = variantInfo[, 1]
    position = as.integer(variantInfo[, 2])
    controlInNonPAR = XYNonPAR(chrID, position, reference)

    variantInfo = matrix(unlist(strsplit(caseCountQC[, 1], "-")), 
      ncol = 4, byrow = T)
    chrID = variantInfo[, 1]
    position = as.integer(variantInfo[, 2])
    caseInNonPAR = XYNonPAR(chrID, position, reference)

    for (i in 1:length(controlAN)) {
      if (extractStratifiedCounts) {
        nCase = ceiling(max(caseACANExtracted[[2]][!caseInNonPAR, i]) / 2)
      } else {
        nCase = ceiling(max(caseCountQC[!caseInNonPAR, caseAN[i]]) / 2)
      }
      nControl = ceiling(max(controlCountQC[!controlInNonPAR, controlAN[i]]) / 2)
      print(paste("group:", ACANConfigData[i, "caseGroup"],
                  "cases:", nCase, "controls:", nControl))
    }
  }
  
  if (!is.null(ACANConfigData)) {
    if (!("caseGroup" %in% colnames(ACANConfigData))) {
      nGroup = 1
    } else {
      nGroup = length(unique(ACANConfigData[, "caseGroup"]))
    }
    # extract AC AN config for controls
    controlACANConfig = ACANConfigData[, c("controlAC", "controlAN"), 
                                       drop = F]
    # make sure it is numeric
    IDs = c(controlACANConfig)
    for (i in 1:length(IDs)) {
      controlCountQC[, IDs[i]] = as.numeric(controlCountQC[, IDs[i]])
    }
    
    # estimate the number of controls per group
    nControlPerGroup = matrix(NA, nGene, nGroup)
    for (i in 1:nGroup) {
      controlANGroup = controlACANConfig[i, 2]
      nControlPerGroup[, i] = ceiling(max(controlCountQC[, controlANGroup], 
                                          na.rm = T) / 2)
    }
    
    if (!is.null(caseGroupInfo)) {
      # leave the extraction of AC AN later using full genotypes
      
      # extract nCase
      nCasePerGroup = matrix(NA, nGene, nGroup)
      for (i in 1:nGroup) {
        nCasePerGroup[, i] = sum(caseGroupInfo == 
                                   ACANConfigData[, "caseGroup"][i])
      }
      
    } else if (("caseAC" %in% colnames(ACANConfigData) && 
                ("caseAN" %in% colnames(ACANConfigData)))) {
      # extract AC AN config for cases
      caseACANConfig = ACANConfigData[, c("caseAC", "caseAN"), drop = F]
      
      # make sure it is numeric
      IDs = c(caseACANConfig)
      for (i in 1:length(IDs)) {
        caseCountQC[, IDs[i]] = as.numeric(caseCountQC[, IDs[i]])
      }
      
      # estimate nCase
      nCasePerGroup = matrix(NA, nGene, nGroup)
      for (i in 1:nGroup) {
        caseANGroup = caseACANConfig[i, 2]
        nCasePerGroup[, i] = ceiling(max(caseCountQC[, caseANGroup]) / 2)
      }
    }
  } else {
    nGroup = 1
    nControlPerGene = matrix(nControl, nGene, 1)
    nCasePerGene = matrix(nCase, nGene, 1)
    nControlPerGroup = matrix(nControl, nGene, 1)
    nCasePerGroup = matrix(nCase, nGene, 1)
    testResult = matrix(NA, length(geneList), 10)
    controlACANConfig = matrix(c(controlAC, controlAN), 1, 2)
    caseACANConfig = matrix(c(caseAC, caseAN), 1, 2)
  }
  
  # assign matched ethnicity for each group when using precomputed LD pairs
  matchedEthnicity = NULL
  if (!ignoreEthnicityInLD) {
    if (gnomADVersion %in% c("v2exome", "v3genome", "v4exome", "v4genome")) {
      matchedEthnicity = character(6)
      ethnicitySetAll = c("nfe", "afr", "amr", "eas", "sas", "fin")
    } else if (gnomADVersion == "v2genome") {
      matchedEthnicity = character(5)
      ethnicitySetAll = c("nfe", "afr", "amr", "eas", "fin")
    } else {
      stop("gnomAD version should be: v2exome, v2genome, v3genome, v4exome, v4genome")
    }
    for (i in 1:dim(controlACANConfig)[1]) {
      for (j in 1:length(ethnicitySetAll)) {
        if (grepl(ethnicitySetAll[j], controlACANConfig[i, 1], fix = T)) {
          matchedEthnicity[i] = ethnicitySetAll[j]
        }
      }
    }
  }
  
  caseGeneList = list()
  controlGeneList = list()
  caseCountEstimated = NULL
  controlCountEstimated = NULL
  
  if (regexpr("^sample", countUse) == -1)  {
    countMatrix = matrix(NA, length(geneList), 4)
  } else {
    # These estimate counts are arranged as 
    # g1_dominant g2_dominant ... g1_recessive ... g1_twohets ...
    caseCountEstimated = matrix(0, length(geneList), nGroup * 3)
    controlCountEstimated = matrix(0, length(geneList), nGroup * 3)
  }
  
  # extract the counts for each gene for the cases
  excludeHighLDInCase = checkHighLDInControl
  if (fullCaseGenotype && excludeHighLDInCase || 
      regexpr("^sample", countUse) != -1) {
    for (i in 1:length(geneList)) {
      gene = geneList[i]
      variantInCase = caseCountQC[
        caseCountQC[, groupColumn] == gene & caseFunction, 1]
      caseGeneList[[i]] = c(gene, variantInCase)
    }
  }
  
  for (i in 1:nGene) {
    gene = geneList[i]
    
    if (countUse %in% c("allele_AN", "allele_null")) {
      # allele_AN is functional allele counts vs rest allele counts
      # allele_null is functional vs nonfunctional, similar as proxecat 
      # This assumes the simple one group case control AC AN configuration
      
      controlCountFunction = sum(controlCountQC[
        controlCountQC[, groupColumn] == gene & controlFunction, controlAC])
      caseCountFunction = sum(caseCountQC[
        caseCountQC[, groupColumn] == gene & caseFunction, caseAC])
      
      controlCountNonFunction = 0
      caseCountNonFunction = 0
      
      # based on AN
      controlANRow = as.numeric(
        controlCountQC[controlCountQC[, groupColumn] == gene, 
                       controlAN])
      if (length(controlANRow) > 0) {
        controlCountNonFunction = round(mean(controlANRow)) - 
          controlCountFunction
        #controlCountNonFunction = sum(controlANRow) - controlCountFunction
        # nControlPerGene[i, ] = ceiling(max(controlANRow) / 2)
      }
      caseANRow = as.numeric(caseCountQC[caseCountQC[, groupColumn] == gene, 
                                         caseAN])
      if (length(caseANRow) > 0) {
        caseCountNonFunction = round(mean(caseANRow)) - caseCountFunction
        #caseCountNonFunction = sum(caseANRow) - caseCountFunction
        # nCasePerGene[i, ] = ceiling(max(caseANRow) / 2)
      }
      
      if (countUse %in% c("allele_AN")) {
        countMatrix[i, ] = c(caseCountFunction, caseCountNonFunction,
                            controlCountFunction, controlCountNonFunction)
      } else if (countUse == "allele_null") {
        # based on synonymous and other defined nonFunctional variants' AC
        controlCountNonFunction = sum(controlCountQC[
        controlCountQC[, groupColumn] == gene & controlNonFunction, controlAC])
        caseCountNonFunction = sum(caseCountQC[
          caseCountQC[, groupColumn] == gene & caseNonFunction, caseAC])
        
        countMatrix[i, ] = c(caseCountFunction, caseCountNonFunction,
                             controlCountFunction, controlCountNonFunction)
      } 
    } else if (regexpr("^sample", countUse) != -1) {
      ## estimate on the controls
      variantInControlSelected = (controlCountQC[, groupColumn] == gene & 
                                    controlFunction)
      variantInControl = controlCountQC[variantInControlSelected, 1]
      controlGeneList[[i]] = c(gene, variantInControl)

      for (j in 1:nGroup) {
        controlACGroup = controlACANConfig[j, 1]         
        controlDominant = sum(controlCountQC[variantInControlSelected,
          controlACGroup])
        variaintControlAF = (controlCountQC[
          variantInControlSelected, 
          controlACGroup]) / (2 * nControlPerGroup[i, j])
        names(variaintControlAF) = variantInControl
        countHomAlt = NULL
        if ("controlAA" %in% colnames(ACANConfigData)) {
          countHomAlt = as.numeric(controlCountQC[
            variantInControlSelected, 
            ACANConfigData[j, "controlAA"] ])
        } 
        
        if (length(variaintControlAF) > 0) {
          highLDList = list(NULL)
          # if (excludeHighLDInControl || verifyHighLDInControl) {
          #   ethnicitySetAll = c("nfe", "afr", "amr", "eas", "sas", "fin")
          #   if (controlACGroup == "AC") {
          #     ethnicitySet = ethnicitySetAll
          #   } else {
          #     ethnicityIndex = grepl(sub("AC_", "", controlACGroup), 
          #                            ethnicitySetAll)
          #     if (sum(ethnicityIndex) == 0) {
          #       stop(paste("cannot identify the ethnicity group from AC ID:",
          #                  controlACGroup))
          #     }
          #     ethnicitySet = ethnicitySetAll[ethnicityIndex]
          #   }
          #   features = c(1, grep("(AC|AN).*(nfe|amr|afr|eas|sas|fin)", 
          #                   colnames(controlCountQC)))
          #   controlDataPerGene = controlCountQC[variantInControlSelected, 
          #                                       features, 
          #                                       drop = F]
          #   highLDList = identifyHighLDVariants(controlDataPerGene, 
          #                                     pThreshold = pLDControl,
          #                                     ORThreshold = ORThresholdControl,
          #                                     ethnicitySet = ethnicitySet,
          #                                     version = gnomADVersion)
          #   if (verifyHighLDInControl) {
          #     controlHighLD[i, j] = highLDList
          #     if (!excludeHighLDInControl) {
          #       highLDList = list(NULL)
          #     }
          #   }
          # } 
          if (!is.null(highLDLoaded)) {
            excludeHighLDInControl = T
            if (ignoreEthnicityInLD) {
              highLDList = list(highLDLoaded[[1]][[gene]])
            } else {
              highLDList = list(
                highLDLoaded[[1]][[gene]][[ matchedEthnicity[j] ]])
            }
          }
          twoHaps = F
          estimatedControlCounts = estimateCounts(nControlPerGroup[i, j], 
                                      variaintControlAF, 
                                      excludeHighLD = excludeHighLDInControl,
                                      highLDList, 
                                      countHomAlt, 
                                      twoHaps = twoHaps)
          controlCountEstimated[i, c(j, nGroup + j, nGroup * 2 + j)] = 
            estimatedControlCounts
        } else {
          controlCountEstimated[i, c(j, nGroup + j, nGroup * 2 + j)] = 
            c(controlDominant, 0, 0) 
        }
      }
    } else {
      stop("wrong countUse")
    }
  }
  
  ## extract high LD variants in cases
  if (fullCaseGenotype && excludeHighLDInCase) {
    caseHighLD = matrix(list(), length(geneList), nGroup)
    
    # extract LD in controls
    # highLDListRelaxed = matrix(list(), length(geneList), nGroup)
    for (i in 1:nGene) {
      gene = geneList[i]
      for (j in 1:nGroup) {
        highLDList = list(NULL)
        if (!is.null(highLDLoaded)) {
          if (ignoreEthnicityInLD) {
            highLDList = list(highLDLoaded[[2]][[gene]])
          } else {
            highLDList = list(
              highLDLoaded[[2]][[gene]][[ matchedEthnicity[j] ]])
          }
        }
        caseHighLD[i, j] = highLDList
      }
    }
    
    # The following could add LD checks in cases, but not with pThreshold = 0
    # indexVariants = which(sapply(caseGeneList, length) > 1)
    # groupIDSpecified = NULL
    # if (!is.null(controlACANConfig) && !is.null(caseGroupInfo)) {
    #   groupIDSpecified = ACANConfigData[, "caseGroup"]
    # }
    # #browser()
    # if (length(indexVariants) > 0) {
    #   if (!checkHighLDInControl) {
    #     # This is not a thorough check in controls, do not use this
    #     caseHighLDExtracted = extractHighLDFromCases(caseGDS, sampleID,
    #                                                caseGeneList[indexVariants],
    #                                                caseGroupInfo,
    #                                                groupIDSpecified,
    #                                                ORThreshold = 1,
    #                                                pThreshold = 0,
    #                                                ACThreshold = 3,
    #                                                NULL)
    #   } else { # use this LD check
    #     if (!is.null(highLDListRelaxed)) {
    #       caseHighLDExtracted = extractHighLDFromCases(caseGDS, sampleID,
    #                              caseGeneList[indexVariants],
    #                              caseGroupInfo,
    #                              groupIDSpecified,
    #                              ORThreshold = 1,
    #                              pThreshold = 0,
    #                              ACThreshold = 3,
    #                              NULL,
    #                              checkHighLDInControl = T,
    #                              gnomADVersion = gnomADVersion,
    #                              controlCountQC = NULL,
    #                              groupColumn = NULL,
    #                              controlFunction = NULL,
    #                              controlACANConfig = NULL,
    #                              pLDControl = pLDControl,
    #                              ORThresholdControl = ORThresholdControl,
    #                              highLDListRelaxed = highLDListRelaxed)
    #     } else {
    #       caseHighLDExtracted = extractHighLDFromCases(caseGDS, sampleID,
    #                            caseGeneList[indexVariants],
    #                            caseGroupInfo,
    #                            groupIDSpecified,
    #                            ORThreshold = 1,
    #                            pThreshold = 0,
    #                            ACThreshold = 3,
    #                            NULL,
    #                            checkHighLDInControl = T,
    #                            gnomADVersion = gnomADVersion,
    #                            controlCountQC = controlCountQC,
    #                            groupColumn = groupColumn,
    #                            controlFunction = controlFunction,
    #                            controlACANConfig = controlACANConfig,
    #                            pLDControl = pLDControl,
    #                            ORThresholdControl = ORThresholdControl)
    #     }
    #   }
    #   caseHighLD[indexVariants, ] = caseHighLDExtracted
    # }
    
  }
  
  if (regexpr("sample", countUse) != -1 && fullCaseGenotype) {
    # extract the counts for each gene for the cases
    caseCountOutput = matrix(0, length(geneList), nGroup * 3)
    indexVariants = which(sapply(caseGeneList, length) > 1)
    groupIDSpecified = NULL
    if (!is.null(controlACANConfig) && !is.null(caseGroupInfo)) {
      groupIDSpecified = ACANConfigData[, "caseGroup"]
    }
    #browser()
    if (length(indexVariants) > 0) {
      caseExtracted = extractCounts(caseGDS, sampleID, 
                                    caseGeneList[indexVariants],
                                    caseGroupInfo, groupIDSpecified,
                                    excludeHighLDInCase, 
                                    caseHighLD[indexVariants, , drop = F])
      caseCountOutput[indexVariants, ] = caseExtracted$extracted
      if (!is.null(caseGroupInfo)) {
        caseCountEstimated[indexVariants, ] = caseExtracted$estimated
      }
    }
  }
  
  # add estimated case counts, this might be useful to check with the 
  # actually extracted counts
  if (regexpr("^sample", countUse) != -1) {
    for (i in 1:nGene) {
      gene = geneList[i]
      if (!is.null(caseACANConfig) && is.null(caseGroupInfo)) {
        variantInCaseSelected = (caseCountQC[, groupColumn] == gene & 
                                   caseFunction)
        variantInCase = caseCountQC[variantInCaseSelected, 1]
        for (j in 1:nGroup) {
          caseACGroup = caseACANConfig[j, 1] 
          caseDominantEstimate = sum(caseCountQC[variantInCaseSelected, 
                                                 caseACGroup])
          variaintCaseAF = (caseCountQC[variantInCaseSelected, 
                                        caseACGroup] + 0.1) / 
            (2 * nCasePerGroup[i, j])
          names(variaintCaseAF) = variantInCase
          countHomAlt = NULL
          if (!is.null(caseAA)) {
            countHomAlt = as.numeric(caseCountQC[variantInCaseSelected, 
                                                 caseAA])
          } 
          highLDList = list(NULL)
          if (!is.null(highLDLoaded)) {
            excludeHighLDInControl = T
            if (ignoreEthnicityInLD) {
              highLDList = list(highLDLoaded[[1]][[gene]])
            } else {
              highLDList = list(
                highLDLoaded[[1]][[gene]][[ matchedEthnicity[j] ]])
            }
          }
          twoHaps = F
          if (length(variaintCaseAF) > 0) {
            estimatedCaseCounts = estimateCounts(nCasePerGroup[i, j], 
                                                 variaintCaseAF,
                                      excludeHighLD = excludeHighLDInControl,
                                                 highLDList, 
                                                 countHomAlt, 
                                                 twoHaps = twoHaps)
            caseCountEstimated[i, c(j, nGroup + j, nGroup * 2 + j)] = 
              c(caseDominantEstimate, estimatedCaseCounts[2:3]) 
          } else {
            caseCountEstimated[i, c(j, nGroup + j, nGroup * 2 + j)] = 
              c(caseDominantEstimate, 0, 0)
          }
        }
      } else if (!is.null(caseGroupInfo)) {
        # extract the AC AN and do the estimation 
        # using the full genotypes in the gds file 
      } else {
        stop("cannot estimate the case counts")
      }
    }
  }
  
  ## compare functional vs nonfunctional allele counts
  if (countUse %in% c("allele_null")) {
    testResult = matrix(NA, length(geneList), 10)
    require(proxecat)
    require(DescTools)
    for (i in 1:length(geneList)) {
      count = rbind(countMatrix[i, 1:2],countMatrix[i, 3:4])
      
      out = proxecat(count[1, 1], count[1, 2], count[2, 1], count[2, 2])
      
      fisherResult = list()
      fisherResult$estimate = NA
      fisherResult$p.value = NA
      try({
        fisherResult = fisher.test(count)
      })
      # chisqTestStat = 0
      # chisqPvalue = 1
      # if (count[1, 1] > 0) {
      #   chisqResult = chisq.test(count, correct = F)
      #   chisqTestStat = chisqResult$statistic
      #   chisqPvalue = chisqResult$p.value
      # }
      GTestStat = NA
      GPvalue = NA
      if (sum(count) > 0) {
        GTestStat = GTest(count)$statistic
        GPvalue = GTest(count)$p.value
      }
      
      testResult[i, ] = c(c(as.matrix(out$dat)), out$test.statistic, 
                          out$p.value, 
                          # chisqTestStat, chisqPvalue,
                          GTestStat, GPvalue,
                          fisherResult$estimate, fisherResult$p.value)
    }

    output = data.frame(geneList, testResult, stringsAsFactors = F)
    colnames(output) = c("gene", "caseFunc", "caseNonFunc", 
                         "controlFunc", "controlNonFunc", 
                         "proxecat.testStat", "proxecat.pvalue", 
                         
                         # "chisq.testStat", "chisq.pvalue",
                         "G.testStat", "G.pvalue",
                         "fisher.OR", "fisher.pvalue"
    )
    output = output[order(output[, "fisher.pvalue"]), ]
  }

  if (regexpr("sample", countUse) != -1) {
    if (!fullCaseGenotype) {
      caseCountOutput = caseCountEstimated
    }
    
    # generate the count table
    groupIDs = ""
    if (nGroup > 1) {
      groupIDs = paste0("_", ACANConfigData[, "caseGroup"])
    }
    dominantTable = generateCountTable(caseCountOutput[, 1:nGroup], 
                                       nCasePerGroup,
                                       controlCountEstimated[, 1:nGroup],
                                       nControlPerGroup, 
                                       groupID = paste0(groupIDs, "_DOM"))
    recessiveTable = generateCountTable(
                               caseCountOutput[, (nGroup + 1):(nGroup * 2)], 
                                        nCasePerGroup,
                         controlCountEstimated[, (nGroup + 1):(nGroup * 2)],
                                        nControlPerGroup, 
                         groupIDs = paste0(groupIDs, "_REC"))
    twoHetsTable = generateCountTable(
                    caseCountOutput[, (nGroup * 2 + 1):(nGroup * 3)], 
                    nCasePerGroup,
                    controlCountEstimated[, (nGroup * 2 + 1):(nGroup * 3)],
                    nControlPerGroup, 
                    groupIDs = paste0(groupIDs, "_2HETS"))
    
    if (nGroup == 1) {
      test = "FET"
    } else {
      test = "CMH"
    }
    
    dominantResult = contingencyTableTest(dominantTable, test)
    colnames(dominantResult) = paste0(colnames(dominantResult), "_DOM")
    recessiveResult = contingencyTableTest(recessiveTable, test)
    colnames(recessiveResult) = paste0(colnames(dominantResult), "_REC")
    twoHetsResult = contingencyTableTest(twoHetsTable, test)
    colnames(twoHetsResult) = paste0(colnames(dominantResult), "_2HET")
    
    
    colnames(caseCountEstimated) = rep("", dim(caseCountEstimated)[2])
    colnames(caseCountEstimated)[1:nGroup] = paste0("caseEstimated", 
                                                    groupIDs, "_DOM") 
    colnames(caseCountEstimated)[(nGroup + 1):(nGroup * 2)] = 
                paste0("caseEstimated", groupIDs, "_REC") 
    colnames(caseCountEstimated)[(2 * nGroup + 1) : (nGroup * 3)] = 
                paste0("caseEstimated", groupIDs, "_2HETS") 
    
    output = data.frame(geneList, 
                        dominantResult, dominantTable,
                        recessiveResult, recessiveTable,
                        twoHetsResult, twoHetsTable, caseCountEstimated)
    output = output[order(output[, "P_DOM"]), ]
    colnames(output)[1] = "gene"
  }  
  
  return(list(testResult = output, 
              caseGeneList = caseGeneList, 
              controlGeneList = controlGeneList))
}

