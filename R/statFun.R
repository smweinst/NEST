# Freedman-Lane versions
statFun.lm = function(X, y, Z = 1, type = "coef", n.cores = 1, seed = NULL, FL = FALSE, n.perm = 999, getNull = TRUE){
  
  # determine what column from summary(lm()) output will need to be extracted to match 'type' (if none of the options below are specified, will default to coef)
  type.ind = ifelse(type=="coef", 1, # summary(lm(X[,v] ~ y + Z))$coefficients["y",1] = coefficient
                    ifelse(type == "tvalue", 3, # summary(lm(X[,v] ~ y + Z))$coefficients["y",3] = t-value
                           ifelse(type == "pvalue", 4, # summary(lm(X[,v] ~ y + Z))$coefficients["y",4] = p-value
                                  1))) # revert to coef if unspecified
  
  # pre-specify the permutations (i.e. do the same permutations of people across vertices)
  if (getNull == TRUE){
    if (!is.null(seed)){ set.seed(seed)}
    perm.ind = lapply(1:n.perm, FUN = function(k){
      sample(1:nrow(X), replace = F)
    })
  }
  
  if (FL == FALSE){ # no freedman-lane (assume independence between y and Z)
    stat.obs = unlist(mclapply(1:ncol(X), FUN = function(v){
      # observed statistic:
      fullmod = lm(X[,v] ~ y + Z)
      obs.stat.v = summary(fullmod)$coefficients["y",type.ind]
      
      return(obs.stat.v)
      
    }, mc.cores = n.cores))
    
    if (getNull == TRUE){ # recursive function used below -- specifying getNull = FALSE so that we don't get a null distribution for the null iterations
      stat.null = lapply(1:n.perm, FUN = function(k){
        statFun.lm(X = X, y = y[perm.ind[[k]]], Z = Z, type = type, n.cores = n.cores, seed = seed, FL = FALSE, getNull = FALSE)
      })
      
      #stat.null = do.call("cbind",null.stat)
      
      return(list(T.obs = stat.obs,
                  T.null = stat.null))
    } else{
      return(stat.obs)
    }

    
  }  else{ # do freedman-lane
    if (getNull == TRUE){
      FL.out = mclapply(1:ncol(X), FUN = function(v){
        # observed statistic:
        fullmod = lm(X[,v] ~ y + Z)
        obs.stat.v = summary(fullmod)$coefficients["y",type.ind]
        
        # under H0:
        reducedmod.v = lm(X[,v] ~ Z)
        reducedmod.resid.v = matrix(resid(reducedmod.v), ncol = 1)
        
        null.stat.v = sapply(1:n.perm, FUN = function(k){
          permX.v.k = reducedmod.resid.v[perm.ind[[k]],] + reducedmod.v$fitted.values # permute residuals + add back fitted values from reduced model
          fullmod.permX.v.k = lm(permX.k ~ y + Z) # refit full model using permuted data
          null.stat.v.k = summary(fullmod.permX.v.k)$coefficients["y",type.ind]
          return(null.stat.v.k)
        })
        
        return(list(obs = obs.stat.v, null = null.stat.v))
        
      }, mc.cores = n.cores)
      
      # observed T(v):
      stat.obs = sapply(FL.out, FUN = function(x){x$obs})
      
      # null T(v) (i.e., maps of null T(v) which will be used to get null distribution of enrichment scores):
      stat.null = do.call("cbind",lapply(FL.out, FUN = function(x){x$null}))
      
      return(list(T.obs = stat.obs,
                  T.null = stat.null))
      
    } else{
      message("specify FL = FALSE if getNull = FALSE")
      break
    }

  }

}


# X = image (n x p)
# y = phenotype (n x 1); should be data.frame with name
# Z = covariates (any dimension); should be data.frame with names
# yZ_df = data.frame containing y (phenotype) and Z (any covariates)
# gam.formula = formula object specifying what should be inputted into gam
# lm.formula = formula object specifying what should be inputted into lm (for getting sign at each location)
# y.in.gam = how phenotype variable will appear in gam model matrix (e.g., "s(age)", "sex", or "s(age):sex")
# y.in.lm = how phenotype will appear in lm model matrix
# X = X
# y = data.frame(age = age)
# Z = data.frame(sex = sex)
# gam.formula = as.formula(X.v ~ Z + s(y, fx = T) + s(y, fx = T, by = Z))
# lm.formula = as.formula(X.v ~ y + Z)
# y.permute: specify names of columns in dat that should be permuted

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
      # permute variable of interest -- need to do this before (not after) getting gam model 
      # xyz.k = xyz
      # xyz.k[,y.in.lm] = xyz.k[perm.ind[[k]],y.in.lm]
       # gam.init.k = gam(formula = gam.formula, data = xyz.k)
      
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
