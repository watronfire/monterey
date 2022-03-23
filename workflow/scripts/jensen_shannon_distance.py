import pandas as pd
import time
from epiweeks import Week
from scipy.spatial.distance import jensenshannon
import numpy as np

def load_metadata( md_loc, pair, verbose=True ):
    if verbose:
        print( "Loading metadata...", end="" )
    starting_time = time.time()
    metadata = pd.read_csv( md_loc, usecols=["accession_id", "date_collected", "pangolin_lineage", "site"], parse_dates=["date_collected"] )
    metadata = metadata.loc[metadata["site"].isin(pair)]
    metadata["week"] = metadata["date_collected"].apply( lambda x: Week.fromdate( x ).startdate() )
    metadata["month"] = metadata["date_collected"].to_numpy().astype('datetime64[M]')
    if verbose:
        print( f"Found {metadata.shape[0]} genomes in {time.time() - starting_time:.1f} seconds" )
    return metadata

def genome_bootstrap( df, func, replicates=100 ):
    result_list = list()
    for _ in range( replicates ):
        trial = df.groupby( "site" ).sample( frac=1, replace=True )
        result_list.append( func( trial ) )
    res = np.quantile( result_list, q=[0.025, 0.5, 0.975] )
    return {k:j for k,j in zip( ["js_5", "js_mean", "js_95"], res )}

def apply_jensenshannon( df ):
    evaluate = df.pivot_table( index="pangolin_lineage", columns="site", values="accession_id", aggfunc="count", fill_value=0 )
    evaluate = evaluate.apply( lambda x: x / x.sum() )
    try:
        return jensenshannon( evaluate.iloc[:,0], evaluate.iloc[:,1], base=2 )
    except IndexError:
        return np.nan

if __name__ == "__main__":
    md = load_metadata( snakemake.input.metadata, snakemake.params.pair_list )
    results = md.groupby( snakemake.params.resolution ).apply( js_bootstrap, func=apply_jensenshannon )
    results["siteA"] = snakemake.params.pair_list[0]
    results["siteB"] = snakemake.params.pair_list[1]
    results.to_csv( snakemake.output.results )
