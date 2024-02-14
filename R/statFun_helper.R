# helper functions for statFun() function

getCmat = function(full.mod.columns, columns.to.test){
  
  Cmat = matrix(0, nrow = length(columns.to.test),
                ncol = length(full.mod.columns),
                dimnames = list(columns.to.test,
                                full.mod.columns))
  diag(Cmat[,columns.to.test]) = 1
  
  return(Cmat)
  
}

# H0: C\beta = 0
## constraint matrix C has dimension q x p, where each row corresponds to one constraint (i.e. one constraint for each beta we are testing to be = 0)

mv_wald_stat = function(columns.to.test, vcov.v, beta.v){
  Cmat = getCmat(full.mod.columns = colnames(vcov.v), columns.to.test = columns.to.test)
  to.invert = Cmat %*% vcov.v %*% t(Cmat)
  middle.term = solve(to.invert)
  Cbeta = Cmat %*% beta.v
  stat = t(Cbeta) %*% middle.term %*% Cbeta
  return(stat)
}

get_vcov = function(X.v, model.matrix.full, MTM.inv, beta.v){
  resid.v = (X.v - model.matrix.full%*%beta.v)
  vcov.out = MTM.inv*sum(resid.v^2)/(nrow(model.matrix.full) - length(beta.v)) # n = number of subjects, p = number of parameters in full model
  return(vcov.out)
}

# gam.mod = gam(X[,1] ~ Z + s(y, fx = T) + s(y, fx = T, by = Z))
# system.time(test.update <- update(object = gam.mod, formula = X[,2] ~ .))
# system.time(test.gam2 <-  gam(X[,2] ~ Z + s(y, fx = T) + s(y, fx = T, by = Z)))
# 
# identical(vcov(test.update), vcov(test.gam2))
# 
# 
# gam.mm = model.matrix(gam.mod)
# 
# # bug in model.matrix: remove duplicated column(s) from model matrix
# gam.mm = gam.mm[,-which(duplicated(t(gam.mm))==T)]
# 
# myvcov = get_vcov(mmat = gam.mm, image.v = X[,1])
# Rvcov = vcov(gam.mod)
# identical(myvcov, Rvcov)