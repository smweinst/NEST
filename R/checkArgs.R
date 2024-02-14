#' checkArgs() helper function to check whether necessary arguments have been inputted in NEST() function
#'@param args list of all parameters required in corresponding statFun() (function for getting location-specific association statistics T(v))
#'@param required.args - character vector with names of required arguments (i.e., what parameters *need* to be named in args list)
#'@param optional.args - character vector with anmes of optional arguments (i.e., what parameters will be set to their default versions if unspecified by user)
#'@export

checkArgs = function(args, required.args, optional.args){

  args.defaulting = optional.args[(which(!optional.args %in% names(args)))]

  args.missing = required.args[which(!required.args %in% names(args))]
  if (length(args.missing) > 0){
    # message("1 or more missing or misspecified args for for statFun function: ")
    for (j in 1:length(args.missing)){
      cat("\n")
      cat(args.missing[j])
      cat("\n")
    }
    return(FALSE)
  }else{
    if (length(args.defaulting) > 0){

      cat("\n")
      cat("Using the following default settings: \n")

      if ("Z" %in% args.defaulting){ # default to no covariates (i.e., intercept-only model w/ Z = 1)
        args$Z = 1
        cat("\n args$Z = 1 \n ---> defaulting to no covariates, just intercept \n")
      }

      if ("type" %in% args.defaulting){ # relevant only for statFun.lm
        args$type = "coef"
        cat("\n args$type = 'coef' \n ---> for statFun = 'lm', defaulting to lm coefficients as T(v)")
      }

      if ("FL" %in% args.defaulting){ # (as of 02/09/2024 this should only be relevant for statFun.lm, not the gam version)
        args$FL = FALSE # no freedman-lane
        cat("\n args$FL = FALSE \n ---> defaulting not to use Freedman-Lane permutation \n")
      }

      if ("getNull" %in% args.defaulting){ # default to get null distribution
        args$getNull = TRUE
        cat("\n args$getNull = TRUE \n ---> defaulting to get a null distribution (i.e. not just getting observed T(v)) \n")
      }

      if ("n.perm" %in% args.defaulting){ # default to get null distribution
        args$n.perm = 999
        cat("\n args$n.perm = 999 \n ---> defaulting to 999 permutations \n")
      }
    }
    return(args)
  }

}
