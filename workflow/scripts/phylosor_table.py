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

    return output_df

if __name__ == "__main__":
    t = load_tree( tree_loc=snakemake.input.tree )
    tl = [i.label for i in t.taxon_namespace]

    md = load_metadata( md_loc=snakemake.input.metadata, tip_labels=tl )

    # TODO: this will be the place to generate between or within-location queries
    name_A = snakemake.params.pair_list[0]
    query_A = md["site"]== name_A
    name_B = snakemake.params.pair_list[1]
    query_B = md["site"]== name_B

    output = phylosor_table( tree=t,
                             metadata=md,
                             queryA=query_A,
                             nameA=name_A,
                             queryB=query_B,
                             nameB=name_B,
                             window=int( snakemake.params.window_size ) )

    output.to_csv( snakemake.output.results )