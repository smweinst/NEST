#' statFun.gam.mvwald() will be called by NEST() if NEST argument statFun=="gam.mvwald"
#'@param X n x p matrix of image measurements (n = number of subjects, p = number of image locations)
#'@param dat matrix including phenotype and any covariates. make sure all columns in dat are names
#'@param gam.formula user-specified formula for gam fit at each vertex (in each gam, outcome will be a different column of X and the right side of the formula should include variables included in dat matrix)
#'@param lm.formula user-specified formula for lm fit at each vertex (purpose is just to get a sign of the multivariate wald statistic obtained from the gam at each vertex, since the wald stat is always positive)
#'@param y.in.gam character specifying how phenotype variable can be identified in gam. if phenotype is included as a linear term in the gam formula, this could just be the variable name. if there's a smooth term for the phenotype, it may be something like s(y)
#'@param y.in.lm character specifying how phenotype variable can be identified in lm. again, might just be the name of the variable, but could be different (e.g., if testing an interaction term)
#'@param y.permute specify names of columns in dat matrix that should be permuted when getNull = TRUE
#'@param n.cores for parallelization, number of cores to use (default is 1 / no parallelization)
#'@param seed optional to set seed
#'@param n.perm number of permutations to conduct for inference. default is 999 (i.e. minimum p-value 1/1000)
#'@param getNull whether to obtain null distribution vs. just get observed map of statistics. default will be TRUE inside NEST function, but the statFun function will then be recursively called to get null distribution and getNull will then switch to FALSE
#'@export
statFun.gam.mvwald = function(X, dat, gam.formula, lm.formula, y.in.gam, y.in.lm, y.permute = NULL, n.cores = 1, seed = NULL, n.perm = 999, getNull = TRUE){

  xyz = cbind(X.v = X[,1], # just use the first image location for building model.matrix (will end up refitting at every location though)
              dat)

  # initial gam (to get model matrix)
  gam.init = gam(formula = gam.formula, data = data.frame(xyz))
  gam.model.matrix = model.matrix(gam.init)

  # weird bug: model matrix is not full rank; remove duplicated column
  if (sum(duplicated(t(gam.model.matrix))) > 0){
    gam.model.matrix = gam.model.matrix[,-which(duplicated(t(gam.model.matrix)))]
  }

  columns.to.test = colnames(gam.model.matrix)[grepl(y.in.gam, colnames(gam.model.matrix), fixed = T)]

  # which columns exist in reduced model.matrix
  # columns.reduced = colnames(model.matrix.gam)[-which(colnames(model.matrix.gam) %in% columns.to.test.gam)]

  # pre-specify the permutations (i.e. do the same permutations of people across vertices)
  if (getNull == TRUE){
    if (is.null(y.permute)){
      message("need to specify which columns of `dat` should be permuted")
      return(NULL)
    }
    if (!is.null(seed)){set.seed(seed)}
    perm.ind = lapply(1:n.perm, FUN = function(k){
      sample(1:nrow(X), replace = F)
    })
  }

  # because multivariate wald statistic is always going to be positive, get sign from an lm
  ## initial lm (to get model that can be easily updated)
  stat.signs.from.lm = unlist(mclapply(1:ncol(X), FUN = function(v){
    xyz[,"X.v"] = X[,v]
    lm.v = lm(lm.formula, data = data.frame(xyz))
    unname(sign(coef(lm.v)[y.in.lm]))
  }, mc.cores = n.cores))

  # matrix multiplication needed (do before entering loop since it'll be the same across all locations):
  # (here, 'M' stands for model matrix, T is transpose)
  MTM = crossprod(gam.model.matrix)
  MTM.inv = solve(MTM)
  MTM.inv.MT = MTM.inv%*%t(gam.model.matrix)
  all.beta = MTM.inv.MT%*%X

  #if (FL == FALSE){ # commenting this part out -- haven't yet implemented version with freedman-lane for gam/multivariate wald
  # observed statistic:
  stat.obs = unlist(mclapply(1:ncol(X), FUN = function(v){
    beta.v = matrix(all.beta[,v], nrow = length(all.beta[,v]))
    vcov.v = get_vcov(X.v = X[,v], model.matrix.full = gam.model.matrix, MTM.inv = MTM.inv, beta.v = beta.v)
    mvwald.stat = mv_wald_stat(columns.to.test = columns.to.test, vcov.v = vcov.v, beta.v = beta.v)
    return(mvwald.stat)
  }, mc.cores = n.cores))
  stat.obs = stat.obs*stat.signs.from.lm

  if (getNull == TRUE){
    stat.null = lapply(1:n.perm, FUN = function(k){

      dat.k = dat
      dat.k[,y.permute] = dat.k[perm.ind[[k]], y.permute]

      statFun.gam.mvwald(X = X,
                         dat = dat.k,
                         # y = y[perm.ind[[k]]], # permute phenotype
                         # Z = Z,
                         gam.formula = gam.formula,
                         lm.formula = lm.formula,
                         y.in.gam = y.in.gam,
                         y.in.lm = y.in.lm, n.cores = n.cores, seed = seed,
                         n.perm = NULL, # no permutations to do with getNull = FALSE
                         getNull = FALSE)

    })
    return(list(T.obs = stat.obs,
                T.null = stat.null))
  } else{
    return(stat.obs)
  }


}
