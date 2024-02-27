import numpy as np

def pval_fun(obs, null_dist):
    K = len(null_dist)
    pval = (1 + np.sum(null_dist > obs)) / (1 + K)
    return pval