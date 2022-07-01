import argparse

import numpy as np
from dendropy import Tree
import pandas as pd
import time

def hill( tree, comA, comB, q=1, rel_then_pool=True ):
    """Calculates the phylogenetic beta diversity using Hill numbers for two communities.
    Parameters
    ----------
    tree : dendropy.Tree
        phylogenetic tree containing taxa from both communities.
    comA : list of str
        set of taxa collected from first community.
    comB : list of str
        set of taxa collected from second community.
    q : float
        Hill number, q = 0 to get species richness, q = 1 (default) to get shannon entropy, q = 2 will give inverse Simpson.
    rel_then_pool : bool
        default is True. Abundance of species are first changed to relative abundance within sites, then pooled into one
        assemblage. If False, sites are pooled first, then change abundance of species to relative abundance.

    Returns
    -------
    q : float
        Hill number used in calculations
    gamma_pd : float
        phylogenetic gamma diversity
    alpha_pd : float
        phylogenetic alpha diversity
    beta_pd : float
        phylogentic beta diversity
    local_similarity : float
        local species overlap. Similar to PhyloSor and bound by [0,1].
    region_similarity : float
        region species overlap. Similar to UniFrac and bound by [0,1]
    """

    df = {
        "name" : [],
        "branch_length" : [],
        "commA" : [],
        "commB" : []
    }

    node_iter = 0
    tree0 = tree.extract_tree()
    for node in tree0.postorder_node_iter():

        if node == tree.seed_node:
            continue

        if node.is_leaf():
            name = node.taxon.label
            node.commA = 1 if name in comA else 0
            node.commB = 1 if name in comB else 0
        else:
            name =  f"node_{node_iter}"
            node_iter += 1
            node.commA = 0
            node.commB = 0
            for child in node.child_node_iter():
                node.commA += child.commA
                node.commB += child.commB
        df["commA"].append( node.commA )
        df["commB"].append( node.commB )
        df["name"].append( name )
        df["branch_length"].append( node.edge_length )
    df = pd.DataFrame( df )
    df = df.loc[(df["commA"]>0)|(df["commB"]>0)]

    if rel_then_pool:
        df["commA"] = df["commA"] / df["commA"].sum()
        df["commB"] = df["commB"] / df["commB"].sum()
    df["total_comm"] = df["commA"] + df["commB"]

    gT = ( df["branch_length"] * df["total_comm"] ).sum()
    community_diversity = df[["commA", "commB"]] / gT

    if q == 1:
        gamma_pd = np.exp( -1 * (df["branch_length"] * (df["total_comm"] / gT) * np.log(df["total_comm"] / gT ) ).sum() )
        alpha_pd = np.exp( -1 * ( community_diversity.multiply( df["branch_length"], axis=0 ) * np.log( community_diversity ) ).sum().sum() ) / 2
        beta_pd = gamma_pd / alpha_pd

        local_similarity = 1 - np.log( beta_pd ) / np.log(2)
        region_similarity = local_similarity
    else:
        exponent = 1 / ( 1 - q )
        gamma_pd = np.power( ( df["branch_length"] * np.power( df["total_comm"] / gT, q ) ).sum(), exponent )
        alpha_a = np.power( community_diversity.loc[community_diversity["commA"]>0,"commA"], q ).multiply( df["branch_length"], axis=0 ).sum()
        alpha_b = np.power( community_diversity.loc[community_diversity["commB"]>0,"commB"], q ).multiply( df["branch_length"], axis=0 ).sum()
        alpha_pd = np.power( alpha_a + alpha_b, exponent ) / 2
        beta_pd = gamma_pd / alpha_pd
        local_similarity = 1 - ( np.power( beta_pd, 1.0 - q ) - 1 ) / ( np.power( 2, 1.0 - q ) - 1 )
        region_similarity = 1 - ( np.power( beta_pd, q - 1.0 ) - 1 ) / ( np.power( 2, q - 1.0 ) - 1 )

    return [
        q,
        gamma_pd,
        alpha_pd,
        beta_pd,
        local_similarity,
        region_similarity
    ]

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
    comA : list
        set of taxa collected from first community.
    comB : list
        set of taxa collected from second community.

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
    blBoth = 0

    tree0 = tree.extract_tree()
    for i in tree0.leaf_nodes():
        if i.taxon.label in comA:
            blA += get_edge_length( i )
            for j in i.ancestor_iter():
                if getattr( j, "comA", False ):
                    break
                elif j.edge.length is not None:
                    j.comA = True
                    blA += j.edge.length
                    if getattr( j, "comB", False ):
                        blBoth += j.edge.length
        elif i.taxon.label in comB:
            blB += get_edge_length( i )
            for j in i.ancestor_iter():
                if getattr( j, "comB", False ):
                    break
                elif j.edge.length is not None:
                    j.comB = True
                    blB += j.edge.length
                    if getattr( j, "comA", False ):
                        blBoth += j.edge.length

    return [blA, blB, blBoth]
    #print( f"blA: {blA}" )
    #print( f"blB: {blB}" )
    #print( f"blBoth: {blBoth}" )
    #print( f"phylosor: {blBoth / (0.5* (blA + blB) ):.2f}" )


def load_tree( tree_loc, verbose=True ):
    """ Loads a newick tree from file, with additional functionality to report how long loading took.

    Parameters
    ----------
    tree_loc : str
        path to newick file.
    verbose : bool
        whether to report how long loading took.

    Returns
    -------
    dendropy.Tree
    """
    if verbose:
        print( "Loading tree...", end="" )
    starting_time = time.time()
    tree = Tree.get( path=tree_loc, schema="newick", preserve_underscores=True )
    if verbose:
        print( f"Done in {time.time() - starting_time:.1f} seconds" )
    return tree


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
    metadata["month"] = metadata["date_collected"].astype( 'datetime64[M]' )
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
    assert all(md_shuffled.pivot_table( index="month", columns=["shuffled"], values="accession_id", aggfunc="count",fill_value=0 ) \
               == md_shuffled.pivot_table( index="month", columns=["site"], values="accession_id", aggfunc="count", fill_value=0 ) ),\
    "shuffle column doesn't contain same number of locations"
    return md_shuffled


def comparison_table( tree, metadata, queryA, nameA, queryB, nameB, window, method, verbose=True ):

    method_func = phylosor if method=="phylosor" else hill

    date_seq = metadata["month"].sort_values().unique()
    if verbose:
        print( f"Performing {len( date_seq )} comparisons." )

    output_df = list()
    for i, month in enumerate( date_seq ):
        start_time = time.time()
        string_rep = np.datetime_as_string( month, unit="D" )
        communityA = set( metadata.loc[(metadata["month"]==month)&queryA, "accession_id"].to_list() )
        communityB = set( metadata.loc[(metadata["month"]==month)&queryB, "accession_id"].to_list() )

        if len( communityA ) + len( communityB ) == 0:
            if verbose:
                print( f"{string_rep} being skipped. No sequences")
            continue

        entry = method_func( tree, communityA, communityB )
        entry.extend( [string_rep, nameA, len( communityA), nameB, len( communityB )] )
        output_df.append( entry )

        if verbose:
            print( f"Completed {i} of {len(date_seq)} comparisons (took {time.time() - start_time:.1f} seconds)." )

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

    t = load_tree( tree_loc=args.tree )
    tl = [i.label for i in t.taxon_namespace]

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

    output = comparison_table( tree=t,
                               metadata=md,
                               queryA=query_A,
                               nameA=name_A,
                               queryB=query_B,
                               nameB=name_B,
                               window=int( args.window_size ),
                               method="hill" if args.hill else "phylosor" )

    output.to_csv( args.output, index=False )
