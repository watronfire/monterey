import argparse

import pandas as pd


def collapse_metadata( md_loc, output ):
    na_countries = ['Antigua and Barbuda', 'Bahamas', 'Barbados', 'Belize',
                    'Bermuda', 'Canada', 'Costa Rica', 'Cuba',
                    'Dominica', 'Dominican Republic', 'El Salvador', 'Guadeloupe', 'Grenada',
                    'USA', 'Guatemala', 'Haiti', 'Honduras',
                    'Jamaica', 'Mexico', 'Panama', 'Saint Barthélemy', 'Saint Kitts and Nevis',
                    'Saint Lucia', 'Saint Martin', 'Saint Vincent and the Grenadines', 'Sint Maarten']
    md = pd.read_csv( md_loc )
    md["site"] = "Other"
    md.loc[md["country"].isin( na_countries),"site"] = md["country"]
    md.loc[md["country"].isin( ["Canada", "Mexico"] ), "site"] = md["division"]
    md.loc[md["country"] == "USA", "site"] = md["location"]
    md.loc[md["division"].isin( ["Guam", "Puerto Rico", "Virgin Islands", "Northern Mariana Islands"] ), "site"] = md["division"]

    early = md.sort_values( "date_collected" ).head(50)["accession_id"].to_list()
    md.loc[md["accession_id"].isin( early ),"site"] = "root"

    md.to_csv( output, index=False )

if __name__ == "__main__":
    parser = argparse.ArgumentParser( description="Adds column to metadata with most-relevant location data" )

    parser.add_argument( "--input", help="location of metadata", required=True )
    parser.add_argument( "--output", help="location to save appended metadata", required=True )

    args = parser.parse_args()

    collapse_metadata( args.input, args.output )