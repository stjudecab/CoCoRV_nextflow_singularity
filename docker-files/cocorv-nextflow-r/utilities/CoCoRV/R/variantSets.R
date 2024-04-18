#' @export
variantSets = function(data, 
                       extraParamJason = "",
                       variantGroupCustom = "", 
                       variantGroup = "annovar_pathogenic",
                       minREVEL = 0.65, 
                       maxAFPopmax = 1) {
  # define the variants of interest in the test
  # defined variant sets are: 
  # "annovar_pathogenic", "annovar_function", "annovar_LOF", "annovar_synonym",
  # "annovar_splicing", "annovar_intron", "annovar_UTR3",
  # "annovar_UTR5", "annovar_UTR"
  # It can also call self defined R functions for variant set, where 
  # variantGroupCustom is the R code file, variantGroup is the R function 
  # defining the set, see below customVariantSets as an example
  
  # input
  # data: a data frame
  # variantGroupCustom: an R code file for self defined variant set
  # variantGroup: either a defined variant set name or a function name for 
  #   self defined variant set in variantGroupCustom
  # minREVEL: minimal REVEL score used in annovar_pathogenic and 
  #    annovar_function
  
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
  if (variantGroup %in% definedSet) {
    if (variantGroup %in% c("annovar_pathogenic", 
                            "annovar_function",
                            "annovar_missense")) {
      # define the set of variants considered as pathogenic 
      
      functional = "ExonicFunc.refGene"
      if (variantGroup == "annovar_missense") {
        functionSet = c("nonsynonymous_SNV")
      } else {
        functionSet = c("frameshift_deletion", 
                        "frameshift_insertion",
                        "stopgain",
                        "nonsynonymous_SNV")
      }
      
      missense = "nonsynonymous_SNV"
      pathogenicScore = "REVEL"
      # eps make it immune to not exact float number representation
      eps = 1e-6
      minREVEL = minREVEL - eps   
      
      outIndex = (data[, functional] %in% functionSet & 
                            ((!(data[, functional] %in% missense)) | 
                               (!is.na(data[, pathogenicScore]) & 
                                  data[, pathogenicScore] >= minREVEL))
      )
      
      if (variantGroup == "annovar_function") { # include the splicing
        outIndex = (outIndex | data[, "Func.refGene"] == "splicing")
      }
    } else {
      functional = "ExonicFunc.refGene"
      if (variantGroup == "annovar_LOF") {
        functionSet = c("frameshift_deletion", "frameshift_insertion", 
                        "stopgain")
      } else if (variantGroup == "annovar_synonym") {
        functionSet = "synonymous_SNV"
      } else if (variantGroup == "annovar_splicing") {
        functionSet = "splicing"
      } else if (variantGroup == "annovar_intron") {
        functional = "Func.refGene"
        functionSet = "intronic,ncRNA_intronic"
      } else if (variantGroup == "annovar_UTR3") {
        functional = "Func.refGene"
        functionSet = "UTR3"
      } else if (variantGroup == "annovar_UTR5") {
        functional = "Func.refGene"
        functionSet = "UTR5"
      } else if (variantGroup == "annovar_UTR") {
        functional = "Func.refGene"
        functionSet = "UTR3,UTR5"
      } else {
        stop("variant group not defined")
      }
      
      outIndex = (data[, functional] %in% functionSet)
      if (variantGroup == "annovar_LOF") {
        outIndex = (outIndex | data[, "Func.refGene"] == "splicing")
      }
    }
  } else { # self defined variant set
    if (variantGroupCustom != "") {
      if (regexpr(".R$|.r$", basename(variantGroupCustom)) != -1) {
        source(variantGroupCustom)
        definedFunction = get(variantGroup)
        if (extraParamJason != "") {
          outIndex = definedFunction(data, extraParamJason)
        } else {
          outIndex = definedFunction(data)
        }
      } else {
        # a table of two columns, tab separated 
        # the first is the variantGroup, the second is the expression
        definedVariantGroups = read.table(variantGroupCustom, header = T, 
          as.is = T, sep = "\t", quote = "")
        expressionString = definedVariantGroups[
                     definedVariantGroups[, 1] == variantGroup, 2]
        if (length(expressionString) == 0) {
          stop("no matched variant group logical string")
        }
        # evaluate the logical string and get the subset
        outIndex = eval(parse(text=as.character(expressionString)),
            envir = data)
      }
    } else {
      stop(paste("self defined variant set file not provided,",
                  "please check the option --variantGroupCustom"))
    } 
  }
  
  if (maxAFPopmax < 1) {
    colID = "AF_popmax"
    # round at the 4th decimal place to be consistent with ANNOVAR values
    AFPopmax = as.numeric(data[, colID])
    outIndex = outIndex & 
      (is.na(AFPopmax) | round(AFPopmax, 4) <= maxAFPopmax + 1e-8)  
  }
  
  return(outIndex)
}

customVariantSets = function(data, extraParamJason = NULL) {
  # an example of a custom function to define variant sets as LOF, 
  # assuming annovar based annotation

  # input
  # data: the data frame where columns are used to define variants of interest
  # extraParamJason: if there are extra parameters or data needed, this 
  #     JASON file provide a way to supply extra parameters.
  
  outIndex = (data[, "ExonicFunc.refGene"] %in% c("frameshift_deletion", 
                                                  "frameshift_insertion",
                                                  "stopgain")) | 
              data[, "Func.refGene"] == "splicing"
  
  return(outIndex)
}

