import argparse
import numpy as np
from dendropy import Tree
import pandas as pd
import time
from subprocess import run
from tempfile import NamedTemporaryFile

def get_edge_length( node ):
    """Gets edge length of node, even if not present.
    Parameters
    ----------
    node : dendropy.Node

    Returns
    -------
    float
        branch length subtending node. Will be 0 for root.
    """
    if node.edge_length is None:
        return 0
    else:
        return node.edge_length

def phylosor( tree, comA, comB ):
    """ Calculates the branch lengths of two communities, as well as the branch lengths that they share.

    Parameters
    ----------
    tree : dendropy.Tree
        phylogenetic tree containing taxa from both communities.
    comA : str
        file containing set of taxa collected from first community.
    comB : str
        file containing set of taxa collected from second community.

    Returns
    -------
    blA : float
        total branch length of all taxa in first community.
    blB : float
        total branch length of all taxa in second community.
    blBoth : float
        total branch length shared by both communities
    """
    blA = 0
    blB = 0

    phylosor_output = run( f"pismo --tree {tree} --commA {comA} --commB {comB}", capture_output=True, text=True, shell=True )
    return [i for i in map(float, phylosor_output.stderr.strip().split( "\n" ))]
    #print( f"blA: {blA}" )
    #print( f"blB: {blB}" )
    #print( f"blBoth: {blBoth}" )
    #print( f"phylosor: {blBoth / (0.5* (blA + blB) ):.2f}" )


def load_metadata( md_loc, tip_labels, verbose=True ):
    """ Loads metadata from md_loc, filters to tips specified by tip_labels, and appends a month columns.

    Parameters
    ----------
    md_loc : str
        path to metadata
    tip_labels : list of str
        metadata will be filtered to these accession IDs.
    verbose: bool
        whether to report how long loading takes.

    Returns
    -------
    pandas.DataFrame
        metadata with `month` column added

    """
    if verbose:
        print( "Loading metadata...", end="" )
    starting_time = time.time()
    metadata = pd.read_csv( md_loc, usecols=["accession_id", "date_collected", "site"], parse_dates=["date_collected"] )
    metadata = metadata.loc[metadata["accession_id"].isin( tip_labels )]
    metadata["month"] = metadata["date_collected"].to_numpy().astype('datetime64[M]')
    if verbose:
        print( f"Done in {time.time() - starting_time:.1f} seconds" )
    return metadata


def generate_date_seq( metadata, query, verbose=True ):
    min_date = metadata.loc[query, "date_collected"].min()
    max_date = metadata.loc[query, "date_collected"].max()
    return pd.date_range( start=min_date, end=max_date )


def add_week_to_metadata( metadata, date_column, week_column ):
    from epiweeks import Week
    return_md = metadata.copy()
    return_md[week_column] = return_md[date_column].apply( lambda x: Week.fromdate( x ).startdate() )
    return return_md


def shuffle_within_location( metadata, location, verbose=True ):
    if verbose:
        print( f"Shuffling tips from {location} for inter-location analysis...", end="" )
    start_time = time.time()
    return_md = add_week_to_metadata( metadata, "date_collected", "week" )

    from numpy.random import choice
    return_md.loc[return_md["site"]==location,"site_shuffled"] = return_md.loc[return_md["site"]==location].groupby( "week" )["site"].transform( lambda x: choice( ["A", "B"], len( x ) ) )
    if verbose:
        print( f"Done in {time.time() - start_time:.1f} seconds" )
    return return_md 


def shuffle_locations( metadata, verbose=True ):
    def shuffle_columns( entry, column ):
        entry["shuffled"] = entry[column].sample( frac=1, replace=False ).to_list()
        return entry
    if verbose:
        print( f"Shuffling locations within tree...", end="" )

    start_time = time.time()
    md_shuffled = metadata.groupby( "month" ).apply( shuffle_columns, column="site" )

    assert md_shuffled.shape[0] == md.shape[0], f"Shuffled dataframe doesn't have the same number of rows (shuffled: {md_shuffled.shape[0]} vs. original: {md.shape[0]})"
    assert not md_shuffled["shuffled"].equals( md_shuffled["site"] ), f"Shuffled column is identical to original column"
    #assert all(md_shuffled.pivot_table( index="month", columns=["shuffled"], values="accession_id", aggfunc="count",fill_value=0 ) \
    #           == md_shuffled.pivot_table( index="month", columns=["site"], values="accession_id", aggfunc="count", fill_value=0 ) ),\
    #"shuffle column doesn't contain same number of locations"
    if verbose:
        print( f"Done in {time.time() - start_time:.1f} seconds" )
    return md_shuffled


def comparison_table( tree, metadata, queryA, nameA, queryB, nameB, window, method, verbose=True ):

    method_func = phylosor

    date_seq = metadata["month"].sort_values().unique()
    if verbose:
        print( f"Performing {len( date_seq )} comparisons." )

    output_df = list()
    for i, month in enumerate( date_seq ):
        start_time = time.time()
        string_rep = np.datetime_as_string( month, unit="D" )
        communityA = set( metadata.loc[(metadata["month"]==month)&queryA, "accession_id"].to_list() )
        communityB = set( metadata.loc[(metadata["month"]==month)&queryB, "accession_id"].to_list() )

        if (len( communityA ) == 0) or (len( communityB ) == 0 ):
            if verbose:
                print( f"{string_rep} being skipped. No sequences")
            continue

        with NamedTemporaryFile( mode="w", delete=False ) as commA_file, NamedTemporaryFile( mode="w", delete=False ) as commB_file:
            commA_file.write( "\n".join( communityA ) )
            commB_file.write( "\n".join( communityB ) )

        entry = method_func( tree, commA_file.name, commB_file.name )
        entry.extend( [string_rep, nameA, len( communityA), nameB, len( communityB )] )
        output_df.append( entry )

        if verbose:
            print( f"Completed {i} of {len(date_seq)} comparisons (took {time.time():.1f} seconds)" )

    headers = {
        "phylosor" : ["blA", "blB", "blBoth", "date", "siteA", "countA", "siteB", "countB"],
        "hill" : ["q" ,"gamma_pd" ,"alpha_pd" ,"beta_pd" ,"local_similarity" ,"region_similarity", "date", "siteA", "countA", "siteB", "countB"]
    }

    output_df = pd.DataFrame( output_df, columns=headers[method] )

    if method == "phylosor":
        output_df["value"] = output_df["blBoth"] / ( 0.5 * (output_df["blA"] + output_df["blB"] ) )
        output_df["value_turn"] = output_df["blBoth"] / output_df[["blA","blB"]].min(axis=1)

    return output_df

if __name__ == "__main__":

    # Command-line argument parsing
    parser = argparse.ArgumentParser( description="Calculates phylosor metric between location pair" )
    parser.add_argument( "--tree", type=str, help="input tree", required=True )
    parser.add_argument( "--metadata", type=str, help="metadata for tips of tree", required=True )
    parser.add_argument( "--pair-list", nargs="+", help="Site pair to compare", required=True )
    parser.add_argument( "--window-size", type=int, default=30, help="Window to calculate phylosor over" )
    parser.add_argument( "--shuffle", action='store_true', help="whether or not to calculate the null model" )
    parser.add_argument( "--output", type=str, help="Location to save output", required=True )
    parser.add_argument( "--hill", action="store_true", help="compare locations using phylogenetic beta diversity using hill numbers (default is to use PhyloSor)" )

    args = parser.parse_args()

    gotree_output = run( f"gotree labels -i {args.tree}", capture_output=True, text=True, shell=True ) 
    tl = gotree_output.stdout.split( "\n" )

    md = load_metadata( md_loc=args.metadata, tip_labels=tl )

    pair_list = args.pair_list

    name_A = pair_list[0]
    name_B = pair_list[1]
    if pair_list[0] == pair_list[1]:
        md = shuffle_within_location( md, pair_list[0] )
        query_A = md["site_shuffled"] == "A"
        query_B = md["site_shuffled"] == "B"
    elif args.shuffle:
        md = shuffle_locations( md )
        query_A = md["shuffled"]==name_A
        query_B = md["shuffled"]==name_B
    else:
        query_A = md["site"]== name_A
        query_B = md["site"]== name_B

    output = comparison_table( tree=args.tree,
                               metadata=md,
                               queryA=query_A,
                               nameA=name_A,
                               queryB=query_B,
                               nameB=name_B,
                               window=int( args.window_size ),
                               method="hill" if args.hill else "phylosor" )

    output.to_csv( args.output, index=False )
