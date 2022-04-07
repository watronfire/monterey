import argparse

from dendropy import Tree
import pandas as pd
import time

def get_edge_length( node ):
    if node.edge_length is None:
        return 0
    else:
        return node.edge_length

def phylosor( tree, comA, comB ):
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
    if verbose:
        print( "Loading tree...", end="" )
    starting_time = time.time()
    tree = Tree.get( path=tree_loc, schema="newick", preserve_underscores=True )
    if verbose:
        print( f"Done in {time.time() - starting_time:.1f} seconds" )
    return tree


def load_metadata( md_loc, tip_labels, verbose=True ):
    if verbose:
        print( "Loading metadata...", end="" )
    starting_time = time.time()
    metadata = pd.read_csv( md_loc, usecols=["accession_id", "date_collected", "site"], parse_dates=["date_collected"] )
    metadata = metadata.loc[metadata["accession_id"].isin( tip_labels )]
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
    return_md = add_week_to_metadata( metadata, "date_collected", "week" )
    md_shuffled = return_md.groupby( "week" ).apply( shuffle_columns, column="site" )

    assert md_shuffled.shape[0] == md.shape[0], f"Shuffled dataframe doesn't have the same number of rows (shuffled: {md_shuffled.shape[0]} vs. original: {md.shape[0]})"
    assert not md_shuffled["shuffled"].equals( md_shuffled["site"] ), f"Shuffled column is identical to original column"
    assert all(md_shuffled.pivot_table( index="week", columns=["shuffled"], values="accession_id", aggfunc="count",fill_value=0 ) \
               == md_shuffled.pivot_table( index="week", columns=["site"], values="accession_id", aggfunc="count", fill_value=0 ) ),\
    "shuffle column doesn't contain same number of locations"
    return md_shuffled


def phylosor_table( tree, metadata, queryA, nameA, queryB, nameB, window, verbose=True ):

    date_seq = generate_date_seq( metadata, queryA|queryB )
    if verbose:
        print( f"Performing {len( date_seq[:-window] )} comparisons with window length {window}" )


    output_df = list()
    for i in range( len( date_seq[:-window] ) ):
        start_time = time.time()
        comparison = (date_seq[i],date_seq[i+30])
        communityA = set( metadata.loc[metadata["date_collected"].between( *comparison )&queryA, "accession_id"].to_list() )
        communityB = set( metadata.loc[metadata["date_collected"].between( *comparison )&queryB, "accession_id"].to_list() )

        if len( communityA ) + len( communityB ) == 0:
            if verbose:
                print( f"{date_seq[i].date} being skipped. No sequences")
            continue

        entry = phylosor( tree, communityA, communityB )
        entry.extend( [date_seq[i].strftime( "%Y-%m-%d"), nameA, len( communityA), nameB, len( communityB )] )
        output_df.append( entry )

        if verbose:
            print( f"Completed {i} of {len(date_seq)} comparisons (took {time.time() - start_time:.1f} seconds)." )

    output_df = pd.DataFrame( output_df, columns=["blA", "blB", "blBoth", "date", "siteA", "countA", "siteB", "countB" ] )
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
    args = parser.parse_args()
    print( args )

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

    output = phylosor_table( tree=t,
                             metadata=md,
                             queryA=query_A,
                             nameA=name_A,
                             queryB=query_B,
                             nameB=name_B,
                             window=int( args.window_size ) )

    output.to_csv( args.output )
