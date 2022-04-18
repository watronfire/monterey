import time
import pandas as pd
from subprocess import run
import argparse
from tempfile import NamedTemporaryFile
import sys

def prune_to_pair( md_loc, id_col, date_col, location_col, tree_loc, pair, output_loc, max_date=None ):
    gotree_output = run( f"gotree labels -i {tree_loc}", capture_output=True, text=True, shell=True )
    tip_labels = gotree_output.stdout.split( "\n" )

    md = pd.read_csv( md_loc, usecols=[id_col,date_col,location_col] )

    # Subset metadata to just sequences in tree
    md = md.loc[md[id_col].isin( tip_labels )]

    # Remove sequences after date if specified
    if max_date:
        md = md.loc[md[date_col]<max_date]

    # Subset metadata to just sequences in pair + Other
    md = md.loc[md[location_col].isin(pair + ["Other"])]

    print( f"Pruning to pair: {', '.join( pair )}" )
    for i in md[location_col].value_counts().iteritems():
        print( f"Keeping {i[1]} tips from {location_col}=={i[0]}" )

    starting_time = time.time()
    print( "Pruning tree...", end="" )
    with NamedTemporaryFile( mode="w" ) as tf:
        tf.write( "\n".join( md[id_col].to_list() ) )
        run( f"gotree prune -r -i {tree_loc} -f {tf.name} -o {output_loc}", shell=True )
        print( f"Done in {time.time()-starting_time:.1f} seconds" )

if __name__ == "__main__":
    parser = argparse.ArgumentParser( description="" )

    # Initialize optional arguments
    parser.add_argument( "-m", "--metadata", required=True )
    parser.add_argument( "-t", "--tree", required=True )
    parser.add_argument( "-i", "--id", required=True )
    parser.add_argument( "-l", "--location-col", required=True )
    parser.add_argument( "-p", "--pair", nargs=2, required=True )
    parser.add_argument( "-o", "--output", required=True )
    parser.add_argument( "--date-col", required=True )
    parser.add_argument( "--max-date" )

    args = parser.parse_args()

    prune_to_pair( md_loc=args.metadata,
                   id_col=args.id,
                   date_col=args.date_col,
                   location_col=args.location_col,
                   tree_loc=args.tree,
                   pair=args.pair,
                   output_loc=args.output,
                   max_date=args.max_date )
