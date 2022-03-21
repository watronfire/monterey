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

def phylosor_table( tree_loc, md_loc, siteA, siteB, window, output_loc ):
    print( "Loading tree...", end="")
    starting_time = time.time()
    tree = Tree.get( path=tree_loc, schema="newick", preserve_underscores=True )
    tip_labels = [i.label for i in tree.taxon_namespace]
    print( f"Done in {time.time() - starting_time:.1f} seconds")

    print( "Loading metadata...", end="" )
    md = pd.read_csv( md_loc, usecols=["accession_id", "date_collected", "site"], parse_dates=["date_collected"] )
    md = md.loc[md["accession_id"].isin( tip_labels )]
    starting_time = time.time()
    print( f"Done in {time.time() - starting_time:.1f} seconds" )

    min_date = md.loc[md["site"].isin([siteA, siteB]), "date_collected"].min()
    max_date = md.loc[md["site"].isin([siteA, siteB]), "date_collected"].max()
    date_seq = pd.date_range( start=min_date, end=max_date )
    print( f"Performing {len( date_seq[:-window] )} comparisons with window length {window}" )

    output_df = list()
    for i in range( len( date_seq[:-window] ) ):
        start_time = time.time()
        comparison = (date_seq[i],date_seq[i+30])
        communityA = set( md.loc[md["date_collected"].between( *comparison )&(md["site"]==siteA), "accession_id"].to_list() )
        communityB = set( md.loc[md["date_collected"].between( *comparison )&(md["site"]==siteB), "accession_id"].to_list() )

        if len( communityA ) + len( communityB ) == 0:
            print( f"{date_seq[i].date} being skipped. No sequences")
            continue

        entry = phylosor( tree, communityA, communityB )
        entry.extend( [date_seq[i].strftime( "%Y-%m-%d"), siteA, len( communityA), siteB, len( communityB )] )
        output_df.append( entry )

        print( f"Completed {i} of {len(date_seq)} comparisons (took {time.time() - start_time:.1f} seconds)." )

    output_df = pd.DataFrame( output_df, columns=["blA", "blB", "blBoth", "date", "siteA", "countA", "siteB", "countB" ] )
    output_df["value"] = output_df["blBoth"] / ( 0.5 * (output_df["blA"] + output_df["blB"] ) )

    output_df.to_csv( output_loc )

if __name__ == "__main__":
    phylosor_table(tree_loc=snakemake.input.tree,
                   md_loc=snakemake.input.metadata,
                   siteA=snakemake.params.pair_list[0],
                   siteB=snakemake.params.pair_list[1],
                   window=int( snakemake.params.window_size ),
                   output_loc=snakemake.output.results )
