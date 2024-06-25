import numpy as np

from checkArgs import check_args
from pvalFun import pval_fun
from enrichScore import enrichScore
from statFun_lm import stat_fun_lm


def NEST(statFun,args, net_maps, one_sided=True, seed=None):
    '''
    Input:
        statFun:
        args: 
        network_maps: 

    Output:
        pval: p-value
        ES: observed enrichment score
        ES_perm: enrichment score for permutations
        running_sum: running sum for observed statistics
    '''
    if not (isinstance(net_maps, list) or (isinstance(net_maps, np.ndarray) and net_maps.ndim == 1)):
        print("net_maps argument should be a list, of 1-d numpy array")
        return None      
    if args['X'].shape[1] != len(net_maps):
        print("The network should be specified as a vector with length equal to number of columns of X.")
        return None
    if not isinstance(net_maps, np.ndarray):
        net_maps = np.array(net_maps)
    is_binary = np.all(np.isin(net_maps, [0, 1]))
    if not is_binary:
        print("net.maps should include only 0's and 1's (1 = in network/ROI; 0 = outside network/ROI)")
        return None
    
    if statFun == "lm":
        required_args = ["X", "y"]
        optional_args = ["Z", "type", "FL", "getNull", "n_perm"]

        args = check_args(args, required_args, optional_args)

        if args:

            statFun_out = stat_fun_lm(X=args['X'], y=args['y'], Z=args['Z'], type=args['type'],
                                     seed=seed, FL=args['FL'], n_perm=args['n_perm'],
                                     getNull=args['getNull'])
        else:
            print('Error. Please check args.')
            return None
    elif statFun == 'gam_mnwald':
        # to be added
        print('No implementation.')
        return None
    
    # For custom methods, the user should at least specify the arguments X, y, and method (method name). It is optional for users to add more arguments to args as needed.
    elif statFun == 'custom':
        required_args = ["X", "y", "method"]
        optional_args = []
        #optional_args = ["Z", "type", "FL", "getNull", "n_perm"]
        args = check_args(args, required_args, optional_args)

        if args:
            statFun_out = args["method"](args)
    else:
        return None

    networks = ['net_'+str(x) for x in range(1)]

    ES_obs,running_sum = enrichScore(statFun_out['T_obs'],net_maps,networks,1,one_sided,False)

    ES_null = {}
    for i in range(args['n_perm']):
        ES_perm, _ = enrichScore(statFun_out['T_null'][i],net_maps,networks,1,one_sided,False)
        for net, val in ES_perm.items():
            if net in ES_null:
                ES_null[net].append(val)
            else:
                ES_null[net]=[val]


    for net in networks:
        pval = pval_fun(ES_obs[net],ES_null[net])


    return pval,ES_obs,ES_null,running_sum