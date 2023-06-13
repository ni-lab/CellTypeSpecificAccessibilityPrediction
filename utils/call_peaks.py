import numpy as np
import pandas as pd
from scipy.stats import poisson

# functions from basenji Github repo
# https://github.com/calico/basenji/blob/30b46b3eb4db413360adca50ccee97658436c86c/bin/basenji_test.py#L365
def ben_hoch(p_values):
    """ Convert the given p-values to q-values using Benjamini-Hochberg FDR. """
    m = len(p_values)

    # attach original indexes to p-values
    p_k = [(p_values[k], k) for k in range(m)]

    # sort by p-value
    p_k.sort()

    # compute q-value and attach original index to front
    k_q = [(p_k[i][1], p_k[i][0] * m // (i + 1)) for i in range(m)]

    # re-sort by original index
    k_q.sort()

    # drop original indexes
    q_values = [k_q[k][1] for k in range(m)]

    return q_values

# https://github.com/calico/basenji/blob/30b46b3eb4db413360adca50ccee97658436c86c/bin/basenji_test.py#L303
def call_peaks(targets, return_full=False):
    # call peaks
    targets_peaks = []
    for ti in range(targets.shape[1]):
        test_targets_dnase_ds_ti = targets[:, ti]
        test_targets_ti_lambda = np.mean(test_targets_dnase_ds_ti, axis=0)
        test_targets_pvals = 1 - poisson.cdf(
          np.round(test_targets_dnase_ds_ti) - 1, mu=test_targets_ti_lambda)
        test_targets_qvals = np.array(ben_hoch(test_targets_pvals))
        test_targets_peaks = test_targets_qvals < 0.01
        targets_peaks.append(test_targets_peaks)

    targets_peaks = np.array(targets_peaks).T
    if return_full:
        return targets_peaks
    else:
        peak_in_num_cell_types = targets_peaks.sum(axis=1)
        return peak_in_num_cell_types