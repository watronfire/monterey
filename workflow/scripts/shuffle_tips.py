import pandas as pd
from dendropy import Tree
from epiweeks import Week
from subprocess import run

def shuffle_tips( tree, metadata, id_col, date_col, map_file, output_loc ):
    t = Tree.get( path=tree, schema="newick", preserve_underscores=True )
    tip_labels = [i.taxon.label for i in t.leaf_node_iter()]
    md = pd.read_csv( metadata, usecols=[id_col, date_col], parse_dates=[date_col] )
    md = md.loc[md[id_col].isin( tip_labels )]
    assert md.shape[0] > 0, f"{id_col} column of metadata doesn't contain any tree tips"

    md["week"] = md[date_col].apply( lambda x: Week.fromdate(x).startdate() )

    def shuffle_columns( entry, column ):
        entry["shuffled"] = entry[column].sample( frac=1 ).to_list()
        return entry

    md_shuffled = md.groupby( "week" ).apply( shuffle_columns, column=id_col )

    assert md_shuffled.shape[0] == md.shape[0], f"Shuffled dataframe doesn't have the same number of rows (shuffled: {md_shuffled.shape[0]} vs. original: {md.shape[0]})"
    assert not md_shuffled["shuffled"].equals( md_shuffled[id_col] ), f"Shuffled column is identical to original column"

    md_shuffled[[id_col, "shuffled"]].to_csv( map_file, index=False, sep="\t", header=False )

    run( f"gotree rename -m {map_file} -i {tree} -o {output_loc}", shell=True)

shuffle_tips( snakemake.input.tree,
              snakemake.input.metadata,
              snakemake.params.id_col,
              snakemake.params.date_col,
              snakemake.output.map_file,
              snakemake.output.shuffled_tree )

