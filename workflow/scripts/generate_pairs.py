import argparse

import pandas as pd
from epiweeks import Week
from matplotlib import pyplot as plt
from matplotlib.ticker import PercentFormatter
from itertools import combinations

def load_metadata( md_loc ):
    usecols = ["accession_id", "date_collected", "country", "division", "site"]
    md = pd.read_csv( md_loc, usecols=usecols, parse_dates=["date_collected"])
    md = md.loc[~md["site"].isna()]
    md = md.loc[~md["site"].isin( ["Other","None"] )]
    md["week"] = md["date_collected"].apply( lambda x: Week.fromdate( x ).startdate() )
    return md


def calculate_summary( md ):
    # sequences per site
    seqs = md["site"].value_counts()

    # fraction of sequences per week per site.
    completeness = md.pivot_table( columns="site", index="week", values="date_collected", aggfunc="count" )
    completeness = completeness.loc[completeness.index > pd.to_datetime( "2020-02-23" )]
    completeness = completeness.count() / completeness.shape[0]

    return_df = pd.concat( [seqs, completeness], axis=1, ignore_index=False )
    return_df.columns = ["sequences", "completeness"]

    return_df["country"] = 'United States'
    return_df.loc[md.loc[md["country"] == "Canada","site"].unique(), "country"] = "Canada"
    return_df.loc[md.loc[md["country"] == "Mexico","site"].unique(), "country"] = "Mexico"
    return return_df


def prepare_graph( summary, min_sequences, min_completeness, graph_loc ):
    fig, ax = plt.subplots( dpi=200, figsize=(6, 4) )
    ax.scatter( "sequences", "completeness", data=summary.loc[summary["country"] == "United States"],
                label="United States" )
    ax.scatter( "sequences", "completeness", data=summary.loc[summary["country"] == "Mexico"],
                label="Mexico" )
    ax.scatter( "sequences", "completeness", data=summary.loc[summary["country"] == "Canada"],
                label="Canada" )
    ax.scatter( "sequences", "completeness", data=summary.loc["San Diego_CA"], label="San Diego" )
    ax.set_ylabel( "Proportion epiweeks sequenced" )
    ax.set_xlabel( "Total sequences" )
    ax.axvline( min_sequences, linestyle="dashed", linewidth=1, color="red" )
    ax.axhline( min_completeness, linestyle="dashed", linewidth=1, color="red" )
    ax.text( 900, 0.95, f">{min_sequences} genomes", ha="right", color="red" )
    ax.text( 1, 0.77, f">{min_completeness:.0%} epiweeks", color="red" )
    ax.set_xscale( "log" )
    ax.yaxis.set_major_formatter( PercentFormatter( 1, 0 ) )
    ax.set_ylim( 0, 1 )
    ax.legend( frameon=False )
    plt.tight_layout()
    plt.savefig( graph_loc )


def generate_pairs( md_loc, min_sequences, min_completeness, output_loc, graph_loc, summary_loc, locations ):
    # Load metadata and add site information.
    md = load_metadata( md_loc )

    # Calculate sequences and completeness per site
    summary = calculate_summary( md )
    summary.to_csv( summary_loc )

    prepare_graph( summary, min_sequences, min_completeness, graph_loc )

    selected = summary.loc[(summary["completeness"]>min_completeness)&(summary["sequences"]>min_sequences)].index
    with open( output_loc, "w" ) as output:
        if locations is not None:
            for loc in locations:
                [output.write( f"{loc},{i}\n" ) for i in selected if i != loc]
        else:
            [output.write( f"{pair[0]},{pair[1]}\n" ) for pair in combinations( selected, 2 )]

if __name__ == "__main__":
    parser = argparse.ArgumentParser( description="Generate location pairs from metadata with reasonable sampling" )

    # Initialize arguments
    parser.add_argument( "--metadata", help="location of metadata", required=True )
    parser.add_argument( "--min-sequences", help="Keep locations with at least this many sequences", type=int, required=True )
    parser.add_argument( "--min-completeness", help="Keep locations with sequences collected from at least this many epiweeks", type=float, required=True )
    parser.add_argument( "--output", help="location to save pairs", required=True )
    parser.add_argument( "--graph", help="location to save diagnostic plot", required=True )
    parser.add_argument( "--summary", help="location to save summary statistics for each location", required=True )
    parser.add_argument( "--locations", nargs="+", help="locations to compare suitable locations to. Leave blank for pairwise", default=None )

    args = parser.parse_args()

    generate_pairs(
        md_loc=args.metadata,
        min_sequences=args.min_sequences,
        min_completeness=args.min_completeness,
        output_loc=args.output,
        graph_loc=args.graph,
        summary_loc=args.summary,
        locations=args.locations
    )
