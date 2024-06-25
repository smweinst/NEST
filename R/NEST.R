#' NEST function
#'@param statFun (to do: add descriptions)
#'@param args .....
#'@param net.maps .....
#'@param one.sided .....
#'@param n.cores .....
#'@param seed .....
#'@param what.to.return .....
#'@export
NEST = function(statFun, args, net.maps, one.sided = TRUE, n.cores = 1, seed = NULL, what.to.return = c("pval")){

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

  if (statFun == "gam.mvwald"){

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
