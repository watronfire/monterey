import argparse
import numpy as np
import pandas as pd
from subprocess import run

def metadata_prune( metadata, tree, kind, sampling, location, missing_loc, output ):

    md = pd.read_csv( metadata, parse_dates=["date_collected"] )

    gotree_output = run( f"gotree labels -i {tree}", capture_output=True, text=True, shell=True )
    tree_tips = gotree_output.stdout.split( "\n" )

    # Identify tips that are in tree but aren't in metadata
    missing = np.setdiff1d( tree_tips, md["accession_id"].to_list() )
    print( f"{len( missing )} of {len( md['strain'] )} tips missing" )


    if kind in ["fraction", "count"]:
        sampling_strat = pd.read_csv( sampling, parse_dates=["month"], index_col=0 )
        sampling_strat["fraction"] = sampling_strat["fraction"].astype( int )
        sampling_strat["count"] = sampling_strat["count"].astype( int )
        sampling_strat = sampling_strat[kind].to_dict()

        tree_md = md.loc[md["accession_id"].isin(tree_tips)&(md["location"]==location)]
        tree_md["month"] = tree_md["date_collected"].astype( "datetime64[M]" )

        sampled = []
        for ts, df in tree_md.groupby( "month" ):
            requested = sampling_strat.get( ts, 0 ) # Not sure if this should be 0 or all sequences.
            requested = min( requested, df.shape[0] )

            if requested == 0:
                continue

            sampled.append( df.sample( n=requested, replace=False ) )
        sampled = pd.concat( sampled )

        missing = np.concatenate((missing, tree_md.loc[~tree_md["accession_id"].isin( sampled["accession_id"] ),"accession_id"].to_list()))

    with open( missing_loc, "w" ) as missing_file:
        missing_file.write( "\n".join( missing ) )

    # Prune tips that are missing
    run( f"gotree prune --tipfile {missing_loc} < {tree} > {output}", shell=True )

if __name__ == "__main__":
    parser = argparse.ArgumentParser( description="" )

    # Initialize optional arguments
    parser.add_argument( "--tree", help="ouput location" )
    parser.add_argument( "--metadata", help="ouput location" )
    parser.add_argument( "--sampling", help="ouput location" )
    parser.add_argument( "--location", help="ouput location" )
    parser.add_argument( "--kind", help="ouput location" )
    parser.add_argument( "--missing", help="ouput location" )
    parser.add_argument( "--output", help="ouput location" )

    args = parser.parse_args()

    metadata_prune(
        args.metadata,
        args.tree,
        args.kind,
        args.sampling,
        args.location,
        args.missing,
        args.output
    )
