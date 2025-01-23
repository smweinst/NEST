#' NEST function
#'@importFrom parallel mclapply
#'@param statFun A string specifying the statistical method to use. Options include:
#'   \describe{
#'     \item{"lm"}{Linear regression model (`statFun.lm`).}
#'     \item{"gam.mvwald"}{Generalized additive model with multivariate Wald tests (`statFun.gam.mvwald`).}
#'     \item{"gam.deltaRsq"}{Generalized additive model for delta adjusted R-squared (`statFun.gam.deltaRsq`).}
#'     \item{"custom"}{Custom statistical function provided via `statFun.custom`.}
#'   }
#'@param args A list of arguments required by the selected `statFun`. Specific arguments depend on the chosen method:
#'   \describe{
#'     \item{"lm"}{Requires `X` (design matrix), `y` (response variable). Optional: `Z` (covariates), `type`, `FL`, `getNull`, `n.perm`.}
#'     \item{"gam.mvwald"}{Requires `X`, `dat`, `gam.formula`, `lm.formula`, `y.in.gam`, `y.in.lm`.}
#'     \item{"gam.deltaRsq"}{Requires `X`, `dat`, `gam.full.formula`, `gam.null.formula`, `lm.formula`, `y.in.gam`, `y.in.lm`.}
#'     \item{"custom"}{Depends on the custom function implementation and should match the arguments specified in `statFun.custom`.}
#'   }
#'@param net.maps A list of vectors representing network maps. Each vector contains 0s (outside network/ROI) and 1s (inside network/ROI). The length of each vector must match the number of columns in `X`.
#'@param one.sided Logical. If `TRUE`, performs one-sided testing. Default is `TRUE`.
#'@param n.cores Integer. The number of cores to use for parallel processing. Default is `1`.
#'@param seed Integer. Random seed for reproducibility. Default is `NULL`.
#'@param what.to.return A character vector specifying the outputs to return. Options include:
#'   \describe{
#'     \item{"pval"}{P-values for each network (default).}
#'     \item{"ES"}{Enrichment scores (observed and null distributions).}
#'     \item{"ES.obs"}{Observed enrichment scores.}
#'     \item{"T.obs"}{Observed test statistics.}
#'     \item{"T.null"}{Null test statistics.}
#'     \item{"everything"}{All of the above outputs.}
#'   }
#'@param statFun.custom A user-defined function for custom statistical analysis. Required if `statFun` is set to `"custom"`. This function should accept arguments passed via `args`.
#'@return A list containing outputs specified in `what.to.return`. Default is p-values for each network.
#'@export
NEST = function(statFun, args, net.maps, one.sided = TRUE, n.cores = 1, seed = NULL, what.to.return = c("pval"), statFun.custom=NULL){

  # check input requirements:
  if (!is.list(net.maps)){
    message("net.maps argument should be a list, of vectors (one for each network)")
    return(NULL)
  }

  if (ncol(args$X) != length(net.maps[[1]])){
    message("each network should be specified as a vector with length equal to number of columns of X.")
    return(NULL)
  }

  if (!identical(as.numeric(sort(unique(unlist(net.maps)))), as.numeric(c(0,1)))){ # check that network maps are all 0s and 1s
    message("net.maps should include only 0's and 1's (1 = in network/ROI; 0 = outside network/ROI)")
    return(NULL)
  }

  # map brain-phenotype associations
  if (statFun == "lm"){
    required.args = c("X", "y")
    optional.args = c("Z", "type", "FL", "getNull","n.perm")

    args = checkArgs(args = args, required.args = required.args, optional.args = optional.args)
    if (!isFALSE(args)){
      statFun.out = statFun.lm(X = args$X,
                               y = args$y,
                               Z = args$Z,
                               type = args$type,
                               n.cores = n.cores,
                               seed = seed,
                               FL = args$FL,
                               n.perm = args$n.perm,
                               getNull = args$getNull)
    } else{
      message("fix args!")
      return(NULL)
    }

  }

  else if (statFun == "gam.mvwald"){

    # check args:
    required.args = c("X","dat","gam.formula","lm.formula","y.in.gam","y.in.lm")
    optional.args = c("n.perm")

    args = checkArgs(args = args, required.args = required.args, optional.args = optional.args)

    if (!isFALSE(args)){
      statFun.out = statFun.gam.mvwald(X = args$X,
                                      dat = args$dat,
                                      gam.formula = args$gam.formula,
                                      lm.formula = args$lm.formula,
                                      y.in.gam = args$y.in.gam,
                                      y.in.lm = args$y.in.lm,
                                      y.permute = args$y.permute,
                                      n.cores = n.cores, seed = seed,
                                      n.perm = args$n.perm,
                                      getNull = TRUE # if doing this inside NEST function, assume testing is being done (if just want to get map, could just use the statFun function directly)
                                      )
    }else{
      message("fix args!")
      return(NULL)
    }

  }

  else if (statFun == "gam.deltaRsq"){
    if (!isFALSE(args)){
      statFun.out = statFun.gam.deltaRsq(X = args$X,
                                      dat = args$dat,
                                      gam.full.formula = args$gam.full.formula,
                                      gam.null.formula = args$gam.null.formula,
                                      lm.formula = args$lm.formula,
                                      y.in.gam = args$y.in.gam,
                                      y.in.lm = args$y.in.lm,
                                      y.permute = args$y.permute,
                                      n.cores = n.cores, seed = seed,
                                      n.perm = args$n.perm,
                                      getNull = TRUE)
    }else{
      message("fix args!")
      return(NULL)
    }
  }

  # Custom function todo: test
  else if (statFun == "custom"){
    if(!is.null(statFun.custom)){
      statFun.out = do.call(statFun.custom, args)
    }else {
       message(("fix statFun.custom"))
       return(NULL)
    }
  }


  ES.list = mclapply(net.maps, FUN = function(net.map){
      ES.obs = enrichScore(stat.map = statFun.out$T.obs,
                  net.map = net.map,
                  one.sided = one.sided,
                  save.detail = F)

      ES.null = lapply(1:args$n.perm, FUN = function(k){
        enrichScore(stat.map = statFun.out$T.null[[k]],
                    net.map = net.map,
                    one.sided = one.sided,
                    save.detail = F)
      })

      return(list(ES.obs = ES.obs,
                  ES.null.dist = ES.null))
  }, mc.cores = n.cores)

  pval.list = mclapply(1:length(net.maps), FUN = function(net){
    pvalFun(obs = ES.list[[net]]$ES.obs, null.dist = ES.list[[net]]$ES.null.dist)
  })

  out = list()

  if ("everything" %in% what.to.return){
    out$pval = pval.list
    out$ES = ES.list
    out$T.obs = statFun.out$T.obs
    out$T.null = statFun.out$T.null
    return(out)
  }else{
    if ("pval" %in% what.to.return){ # default is to just return p-values
      out$pval = pval.list
    }

    if ("ES" %in% what.to.return){
      out$ES = ES.list
    }

    if ("ES.obs" %in% what.to.return){
      out$ES.obs = sapply(ES.list, FUN = function(x){x$ES.obs})
    }

    if ("T.obs" %in% what.to.return){
      out$T.obs = statFun.out$T.obs
    }

    if ("T.null" %in% what.to.return){
      out$T.null = statFun.out$T.null
    }
    return(out)
  }


}
