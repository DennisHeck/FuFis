import pandas as pd
import numpy as np
from scipy.stats import gmean


def mor_df(df):
    """
    Takes a df with raw counts and normalises them via median of ratios. Assumes that the columns are samples.
    First, a pseudo-reference is created by taking the row-wise geometric mean. To get the scaling factor per sample,
    we take the median of the ratios of the sample's entries versus the pseudo-reference. Zeroes are removed for the
    geometric mean, and rows with all zeroes will be erased from the df.
    """
    # Remove rows where everything is 0.
    df = df.loc[(df != 0).any(axis=1)]
    # We skip 0s for the geometric mean.
    pseudo_reference = {k: gmean([x for x in df.loc[k].to_list() if x > 0]) for k in df.index}
    median_ratios = {c: np.median([val / pseudo_reference[g] for g, val in df[c].iteritems()]) for c in df.columns}
    nonzero_ratios = [k for k, ratio in median_ratios.items() if ratio > 0]
    df = df[nonzero_ratios]
    norm_df = df.apply(lambda x: x / median_ratios[x.name], axis=0)
    return norm_df
