import numpy as np
import statsmodels.api as sm
from statsmodels.gam.api import GLMGam, BSplines

def statFun_gam_deltaRsq(X, dat, gam_full_smoother, gam_null_smoother, y_in_gam, y_in_lm, y_permute=None, seed=None, n_perm=999, get_null=True):
    if seed is not None:
        np.random.seed(seed)
    print('there')
    # Function to calculate delta adjusted R-squared for each column
    def calculate_delta_rsq(v):
        print(v)
        X_v = X[:, v]
        xyz = np.column_stack([X_v, dat])

        # Fit the full model (with smooth term)
        gam_fullmodel = GLMGam(X_v, smoother=gam_full_smoother, exog=dat[[y_in_gam]]).fit()
        print(v,'.1')
        # Fit the null model (without smooth term)
        if gam_null_smoother == None:
            gam_nullmodel = sm.GLM(X_v, dat[[y_in_gam]]).fit()
        else:
            gam_nullmodel = GLMGam(X_v, smoother=gam_null_smoother, exog=dat[[y_in_gam]]).fit()
        print(v,'.2')
        # Calculate delta adjusted R-squared
        adj_rsq = gam_fullmodel.pseudo_rsquared() - gam_nullmodel.pseudo_rsquared()

        # Fit linear model to determine the sign using the variable specified in `y_in_lm`
        lm_model = sm.OLS(X_v, dat).fit()

        print(v,'.3')
        sign_t = np.sign(lm_model.params[dat.columns.get_loc(y_in_lm)])  # Coefficient for the phenotype variable

        # Adjust the sign of delta R-squared
        adj_rsq *= sign_t

        return adj_rsq

    # Compute for each column in X
    stat_obs = np.array([calculate_delta_rsq(v) for v in range(X.shape[1])])

    # If get_null is True, calculate the null distribution
    if get_null:
        permuted_stats = []
        for k in range(n_perm):
            perm_ind = np.random.permutation(X.shape[0])
            dat_k = dat.copy()
            dat_k[y_permute] = dat_k.iloc[perm_ind, dat_k.columns.get_loc(y_permute)]

            # Recompute the statistic for the permuted data
            permuted_stat = statFun_gam_deltaRsq(X, dat_k, gam_full_smoother, gam_null_smoother, y_in_gam, y_in_lm, y_permute=None, seed=seed, get_null=False)
            permuted_stats.append(permuted_stat)

        return {"T_obs": stat_obs, "T_null": permuted_stats}

    return stat_obs