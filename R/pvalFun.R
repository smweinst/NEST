#' pvalFun() helper function to calculate a p-value based on user-specified observed statistic and null distribution
#'@param obs - observed statistic (i.e. enrichment score)
#'@param null.dist - null distribution (i.e., enrichment scores under H0, obtained via permutation)
#'@export

pvalFun = function(obs, null.dist){
  K = length(null.dist)
  pval = (1+length(which(null.dist > obs)))/(1+K)
  return(pval)
}
