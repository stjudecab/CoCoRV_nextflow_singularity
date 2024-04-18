#' Computing p-values and their support for Cochran-Mantel-Haenszel (CMH) 
#' exact tests
#'
#' Computes discrete raw p-values and their support
#' for Cochran-Mantel-Haenszel (CMH) exact tests when applied to 
#' 2 x 2 x K contingency tables.
#'
#' Assume that we have K stratified 2x2 contingency tables resulting from two
#' random variables, both have only two values. The CMH assumes a common 
#' odds ratio for all K 2x2 contingency tables. Under the null hypothesis, 
#' the common odds ratio is 1, i.e., no association between the two random 
#' variables after stratified by another variable which has K values. 
#'
#' To calculate the probability of support, the internal R C function C_d2x2xk
#' used by mantelhaen.test is also used here.
#'
#' The code for the computation of the p-values of CMH
#' exact test is inspired by the source code of mantelhaen.test and the function
#' \code{\link[DiscreteFDR]{fisher.pvalues.support}}.
#'
#' See the Wikipedia article about Fisher's
#' exact test, paragraph Example, for a good depiction
#' of what the code does for each possible value
#' of \code{alternative}.

#' @param counts        a matrix of 4xK columns and any number of lines,
#'                     each line representing a 2 x 2 x K contingency table to
#'                     test. Make sure all 2x2 contingency tables are stored 
#'                     in the same order into a 4 element vector, which are 
#'                     concatenated to form a row with length of 4xK. 
#' @param alternative   same argument as in \code{\link{fisher.test}}. The three
#'                     possible values are \code{"greater"} (default),
#'                     \code{"two.sided"} or \code{"less"} and you can specify
#'                     just the initial letter.

#' @return
#' A list of two elements:
#' \item{raw}{raw discrete p-values.}
#' \item{support}{a list of the supports of the CDFs of the p-values.
#' Each support is represented by a vector in increasing order.}

#' @export
cmh.pvalues.support = function(counts, alternative = "two.sided"){
  stopifnot(alternative %in% c("two.sided", "less", "greater"))
  
  nRow = nrow(counts)
  nColumn = ncol(counts)
  
  pCDFlist = vector("list", nRow)
  raw.pvalues = numeric(nRow)
  
  if (nColumn %% 4 != 0 || nColumn == 0) {
    stop("The number of column of input matrix must be 4xK")
  }
  K = nColumn / 4
  
  for (i in 1:nRow) {
    x = array(counts[i, ], dim = c(2, 2, K))
    
    # Similar as used in mantelhaen.test
    # m: the first column sum, n: the second column sum, of all K tables
    # t: the first row sum of all K tables 
    # s: the sum of the left top element of all K tables
    mn = apply(x, c(2L, 3L), sum)
    m = mn[1L, ]
    n = mn[2L, ]
    t = apply(x, c(1L, 3L), sum)[1L, ]
    s = sum(x[1L, 1L, ])
    lo = sum(pmax(0, t - n))
    hi = sum(pmin(m, t))

    support = lo : hi
    # index of the observed data
    index = which(support == s)
    
    # Density of the distribution on its support
    # The sum of the left top element is used to rank all different 2x2xK
    # table combinations.
    dc = .Call(stats:::C_d2x2xk, K, m, n, t, hi - lo + 1L, PACKAGE = "stats")
    
    if (alternative == "greater") {
       # work on the reversed probability, and keep the reversed 
       # because we want to have pCDFlist[[i]] in increasing order
       # pmin/pmax below is to account for machine rounding issues
       reverved = cumsum(rev(dc))
       pCDFlist[[i]] = pmax(0, pmin(1, reverved))
       reversedIndex = length(support) + 1 - index
       raw.pvalues[i] = pCDFlist[[i]][reversedIndex]
    } else if (alternative == "less") {
       # pmin/pmax below is to account for machine rounding issues
       pCDFlist[[i]] = pmax(0, pmin(1, cumsum(dc)))
       raw.pvalues[i] = pCDFlist[[i]][index]
    } else {
       # ensure that probabilities sum up to 1 
       # (sometimes, needs to be done multiple times)
       atoms= dc
       newsum = sum(atoms)
       while(newsum < 1){
         oldsum = newsum
         atoms = atoms/newsum
         newsum = sum(atoms)
         if(oldsum == newsum) break;
       }
       # pmin/pmax below is to account for machine rounding issues
       # relErr: not to miss equal probabilities due to machine precision
       # relErr = 1 + 1e-7 
       
       # this will be slow when the counts are large
       # pCDFlist[[i]] = pmax(0, pmin(1, 
       #                        sapply(1:length(support), function(nu){
       #                          sum(atoms[atoms <= atoms[nu] * relErr])
       #                          })))
       # raw.pvalues[i] = pCDFlist[[i]][index]
       # # we want to have pCDFlist[[i]] in increasing order:
       # pCDFlist[[i]] = sort(pCDFlist[[i]])
       
       # faster version
       # raw.pvalues[i] = sum(atoms[ atoms <= atoms[index] ])
       # atoms = sort(atoms)
       indexSorted = order(atoms)
       atoms = atoms[indexSorted]
       twoSidedPvalues = twoSidedP(atoms)
       pCDFlist[[i]] = pmax(0, pmin(1, twoSidedPvalues))
       raw.pvalues[i] = pCDFlist[[i]][indexSorted == index]
    }
  }
 
  return(list(raw = raw.pvalues, support = pCDFlist))
}
