import pandas as pd
import os
import argparse


STRUCTURE = "results/phylosor/{kind}.{num}/{pair}/{pair}.actual.csv"
def combine_results( results, output ):
    line_count = 0
    output_df = []
    for result in results:
        result_df = pd.read_csv( result )

        name_split = result.split( "/" )
        assert len( name_split ) == 5, f"{result} is oddly formatted."
        kind_num = name_split[2].split( "." )
        result_df["kind"] = kind_num[0]
        result_df["num"] = kind_num[1]
        
        line_count += result_df.shape[0]
        
        output_df.append( result_df )

    output_df = pd.concat( output_df )

    assert output_df.shape == output_df.drop_duplicates().shape, "columns are duplicated. Need to make more unique"
    assert output_df.shape[0] == line_count, f"Entries dropped from concattenated file, {line_count} counted vs. {output_df.shape[0]} present"

    output_df.to_csv( output, index=False )

if __name__ == "__main__":
    parser = argparse.ArgumentParser( description="Combines all csvs specified in input while appending a column with a summary of the original file name" )
    parser.add_argument( "input", nargs="+", help="input csvs seperated by whitespace" )
    parser.add_argument( "output", help="location to save concatenated file" )
    args = parser.parse_args()

    combine_results( results=args.input, output=args.output )
