import gzip
import itertools
import pandas as pd
from multiprocessing import Pool
from timeit import default_timer as clock

"""Custom function for reading large matrices where pandas has trouble with."""


def split_line(row, separator, string_comparison, string_comparison_start):
    if string_comparison:
        return [x == string_comparison if i >= string_comparison_start else x
                for i, x in enumerate(row.strip().split(separator))]
    else:
        return row.strip().split(separator)


def file_read(file_path, sep='\t', header=0, cores=1, batch_rows=100000, stop_row=None, df_dtypes=None,
              string_comparison=None, string_comparison_start=None):
    """
    Read in a file and process batches of rows in parallel. Afterwards converts everything to a pandas DataFrame.
    Can work better if there are many columns in a file, which pandas doesn't handle well. The file can be non-compressed or gzipped.

    Args:
        sep: Separator where the rows will be split.
        header: Index of the row that will be used as header, every row before will be ignored. Set to None to not have use a header.
        cores: Number of cores that will work on splitting the rows.
        batch_rows: How many rows will be in one batch.
        stop_row: Row index where the read-in will be stopped.
        df_dtypes: Dictionary with {column name: dtype} to enforce certain datatypes for the columns.
        string_comparison: Very specific use case, if your data is binary with '0' and '1' or you want to make it binary by comparing to a specific string. This is faster than
            doing the conversion with df_types. Set string_comparison to the string the entries should be compared to, e.g. '1'.
        string_comparison_start: Index from which on the string comparison will be done, e.g. 4 to have all entries[4:] being compared to the string.

    Returns:
        - Pandas DataFrame of the read table.
    """
    start = clock()
    batch_collector = []
    batch_dfs = []
    header_cols = None
    bin_start = clock()

    def batch_handler(batch_collector):
        """To avoid redundance when processing the last rows that didn't fill a batch."""
        process_pool = Pool(processes=cores)
        row_collector = process_pool.starmap(split_line, zip(batch_collector, itertools.repeat(sep),
                                                             itertools.repeat(string_comparison),
                                                             itertools.repeat(string_comparison_start)))
        batch_df = pd.DataFrame(row_collector, columns=header_cols)
        if df_dtypes:
            batch_df = batch_df.astype(df_dtypes)
        batch_dfs.append(batch_df)

    with gzip.open(file_path, 'rt') if file_path.endswith('.gz') else open(file_path) as file_in:
        for e, entry in enumerate(file_in):
            if e == header:
                header_cols = entry.strip().split('\t')
            elif header is not None and e > header:
                batch_collector.append(entry)
            if e > 0 and len(batch_collector) % batch_rows == 0:
                batch_handler(batch_collector)
                print('batch', e, clock() - bin_start)
                batch_collector = []
                bin_start = clock()
            if stop_row and e > stop_row:
                break

    if len(batch_collector) > 0:  # Wasn't emptied then.
        batch_handler(batch_collector)

    df = pd.concat(batch_dfs)  # In small tests concatenating in the end was faster for some reason.
    print(clock() - start)
    return df

