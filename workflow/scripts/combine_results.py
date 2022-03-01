import pandas as pd
import os

def combine_results( results, output ):
    line_count = 0
    output_df = []
    for result in results:
        result_df = pd.read_csv( result, parse_dates=["date"], index_col=0 )

        name = os.path.splitext( os.path.split( result )[1] )[1]
        name_split = name.split( "." )
        result_df["kind"] = name_split[1]
        result_df["num"] = name_split[2]
        output_df.append( result_df )

    output_df = pd.concat( output_df )

    assert output_df.shape == output_df.drop_duplicates().shape, "columns are duplicated. Need to make more unique"

    output_df.to_csv( output, index=False )

if __name__ == "__main__":
    results = snakemake.input.results_nulls
    results.extend( snakemake.input.results_actual )
    combine_results( results=results, output=snakemake.output.results )