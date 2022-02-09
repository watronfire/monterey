import time
import pandas as pd
from subprocess import run
import argparse
from tempfile import NamedTemporaryFile
import sys

def prune_to_pair( md_loc, id_col, location_col, tree_loc, pair, output_loc ):
    gotree_output = run( f"gotree labels -i {tree_loc}", capture_output=True, text=True, shell=True )
    tip_labels = gotree_output.stdout.split( "\n" )

    md = pd.read_csv( md_loc, usecols=[id_col,location_col] )

    # Subset metadata to just sequences in tree
    md = md.loc[md[id_col].isin( tip_labels )]

    # Subset metadata to just sequences in pair + Other
    md = md.loc[md[location_col].isin(pair + ["Other"])]

    print( f"Pruning to pair: {', '.join( pair )}" )
    for i in md[location_col].value_counts().iteritems():
        print( f"Keeping {i[1]} tips from {location_col}=={i[0]}" )

    starting_time = time.time()
    print( "Pruning tree...", end="" )
    with NamedTemporaryFile( mode="w" ) as tf:
        tf.write( "\n".join( md[id_col].to_list() ) )
        run( f"gotree prune -r -i {tree_loc} -f {tf.name} -o {output_loc}" )
        print( f"Done in {time.time()-starting_time:.1f} seconds" )

with open( snakemake.log[0], "w" ) as f:
    sys.stderr = sys.stdout = f
    prune_to_pair( md_loc=snakemake.input.metadata,
                   id_col=snakemake.params.id_col,
                   location_col=snakemake.params.location_col,
                   tree_loc=snakemake.input.tree,
                   pair=snakemake.params.pair_list,
                   output_loc=snakemake.output.pruned_tree )

#if __name__ == "__main__":
#    parser = argparse.ArgumentParser( description="" )
#
#    # Initialize optional arguments
#    parser.add_argument( "-m", "--metadata", help="ouput location" )
#    parser.add_argument( "-t", "--tree", help="modifier" )
#    parser.add_argument( "-i", "--id", help="modifier" )
#    parser.add_argument( "-l", "--location", help="modifier" )
#    parser.add_argument( "-p", "--pair", nargs=2, help="modifier" )
#    parser.add_argument( "-o", "--output", help="modifier" )
#
#    args = parser.parse_args()
#
#    prune_to_pair( md_loc=args.metadata,
#                   id_col=args.id,
#                   location_col=args.location,
#                   tree_loc=args.tree,
#                   pair=args.pair,
#                   output_loc=args.output)