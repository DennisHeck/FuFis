from collections import Counter
import pandas as pd


def df_column_binner(df, col, num_bins, string_precision=0, lower_bound=None, upper_bound=None, tag='', numerate=True):
    """Takes a pandas DataFrame and returns it with a column with a string indicating to which bin the
    chosen column belongs to, including the size of that bin. When specifying lower and upper_bound, everything
    beyond those bounds will be put into a general group of < / > bound, and only the range between bounds will be
    binned. The tag is added as prefix to new column name."""
    to_bin = df[col].values.tolist()
    if not lower_bound:
        lower_bound = min(to_bin)
    if not upper_bound:
        upper_bound = max(to_bin)
    bin_size = (upper_bound - lower_bound) / num_bins

    bin_idx = []
    for entry in to_bin:
        if entry < lower_bound:
            bin_idx.append([entry, 'lower'])
        elif entry >= upper_bound:
            bin_idx.append([entry, 'upper'])
        else:
            bin_idx.append([entry, (entry-lower_bound)//bin_size])

    bin_occs = Counter([x[1] for x in bin_idx])
    bin_strings = {k: str(round(k*bin_size+lower_bound, string_precision)) + '≤ ' + col + ' <' +
                      str(round((k+1)*bin_size+lower_bound, string_precision)) + (' (#' + str(val) + ')') * numerate
                   for k, val in bin_occs.items() if k not in ['lower', 'upper']}
    bin_strings['lower'] = str(round(lower_bound, string_precision)) + '> ' + col + (' (#' + str(len([x for x in to_bin if x < lower_bound])) + ')')*numerate
    bin_strings['upper'] = str(round(upper_bound, string_precision)) + '≤ ' + col + (' (#' + str(len([x for x in to_bin if x >= upper_bound])) + ')')*numerate

    bin_out = [[x[0], bin_strings[x[1]]] for x in bin_idx]
    df[tag + 'binned ' + col] = [x[1] for x in bin_out]
    return df


def get_distance_to_one(start, end, other):
    """ Assumes that start - end is one region, and gives you the distance to other. other is supposed to be
     1-based, usually a TSS from a TSS file. Returns 0 if other is within the region."""
    distance = min([abs(other - 1 - start), abs(other - end)])  # For start subtract 1 to match 0-based.
    if start < other < end or end < other < start:
        distance = 0
    return distance


