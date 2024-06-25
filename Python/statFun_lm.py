import numpy as np
from sklearn.linear_model import LinearRegression


def stat_fun_lm(X, y, Z=1, type='coef', seed=None, FL=False, n_perm=999, getNull=True):


    y = y.reshape(-1,1)
    Z = Z.reshape(-1,1)
    n, p = X.shape

    if getNull:
        if seed is not None:
            np.random.seed(seed)
        perm_indices = [np.random.permutation(n) for _ in range(n_perm)]

    if not FL:
        # No Freedman-Lane: assume independence between y and Z


        # Combine y and Z for the full model
        # do X = beta* y + gamma * Z regression
        y_stack = np.hstack([y, Z])
        
        # Step 1: Fit the full model
        fullmod = LinearRegression().fit(y_stack, X)
        stat_obs = fullmod.coef_[:,0]
        
        
        if (getNull):
            stat_null = []
            for i in range(n_perm):
                y_perm = y[perm_indices[i]]
                y_stack = np.hstack([y_perm, Z])
                fullmod_perm = LinearRegression().fit(y_stack, X)
                stat_null.append(fullmod_perm.coef_[:,0])
            stat_null = np.array(stat_null)

            return {'T_obs':stat_obs,'T_null':stat_null}
        else:
            return {'T_obs':stat_obs}
        
    else:  # do freedman-lane
        y_stack = np.hstack([y, Z])
        

        y_stack = np.hstack([y, Z])
        
        # Step 1: Fit the full model
        fullmod = LinearRegression().fit(y_stack, X)
        # Compute statistic of interest: coefficient of the variable in y
        stat_obs = fullmod.coef_[:,0]
        
        
        # Step 2: Fit the reduced model (only Z)
        reg_reduced = LinearRegression().fit(Z, X)
        # Compute residuals from the reduced model
        residuals = X - reg_reduced.predict(Z)
        
        stat_null = []
        
        for i in range(n_perm):        
            # Step 3: Generate permuted X* by permuting residuals and adding back the fitted values from Z
            X_perm = reg_reduced.predict(Z) + residuals[perm_indices[i],:]
            
            # Step 4: Regress the permuted X* on the full model
            fullmod_perm = LinearRegression().fit(y_stack, X_perm)
            stat = fullmod_perm.coef_[:,0]  # Statistic from permuted data: null T(v)
            
            stat_null.append(stat)

        stat_null = np.array(stat_null)
        return {'T_obs':stat_obs,'T_null':stat_null}
