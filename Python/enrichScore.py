import numpy as np

def enrichScore(L, L_network_labels,  networks=['test'], p=1, save_vertex_level=False):
    """
    Calculate the enrichment score for each network.

    Parameters:
    L (array-like): Unsorted statistics at each vertex/location.
    L_network_labels (array-like): Labels of networks corresponding to each location. (Binary so far)
    p (int): Exponent in enrichment score calculation. (1)
    networks (list): Network labels for which enrichment score is calculated. (not used)
    save_vertex_level (bool): Whether to save detailed output at the vertex level.

    Returns:
    dict: Enrichment scores for each network.
    """

    # Order locations by statistic
    L_order = np.argsort(L)[::-1]
    L_sorted = L[L_order]

    L_sorted_p = np.abs(L_sorted) ** p
    L_network_labels_sorted = L_network_labels[L_order]
    V = len(L)

    out = {}
    # We can assume there's only one element (1) in networks since we are using binary labels.
    for val, curr_net in enumerate(networks):
        # val is from 0,1,2,..., to len(networks)-1

        curr_net_ind = (L_network_labels_sorted == val+1)

        P_hit_numerator = np.cumsum(L_sorted_p * curr_net_ind)
        P_hit = P_hit_numerator / P_hit_numerator[-1]

        P_miss_numerator = np.cumsum(1 - curr_net_ind)
        P_miss = P_miss_numerator / P_miss_numerator[-1]

        running_sum = P_hit - P_miss 
        ES = np.max(np.abs(running_sum))

        # I comment this part since it's optional and has not debugged yet.
        # if save_vertex_level:
        #     out[curr_net] = {
        #         'P_hit': P_hit,
        #         'P_miss': P_miss,
        #         'net_ind': curr_net_ind,
        #         'ES': ES,
        #         'running_sum': running_sum
        #     }
        # else:
        #     out[curr_net] = ES
        out[curr_net] = ES

    return out,running_sum
