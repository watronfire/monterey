import argparse

import pandas as pd

us_state_to_abbrev = {
    "Alabama": "AL",
    "Alaska": "AK",
    "Arizona": "AZ",
    "Arkansas": "AR",
    "California": "CA",
    "Colorado": "CO",
    "Connecticut": "CT",
    "Delaware": "DE",
    "Florida": "FL",
    "Georgia": "GA",
    "Hawaii": "HI",
    "Idaho": "ID",
    "Illinois": "IL",
    "Indiana": "IN",
    "Iowa": "IA",
    "Kansas": "KS",
    "Kentucky": "KY",
    "Louisiana": "LA",
    "Maine": "ME",
    "Maryland": "MD",
    "Massachusetts": "MA",
    "Michigan": "MI",
    "Minnesota": "MN",
    "Mississippi": "MS",
    "Missouri": "MO",
    "Montana": "MT",
    "Nebraska": "NE",
    "Nevada": "NV",
    "New Hampshire": "NH",
    "New Jersey": "NJ",
    "New Mexico": "NM",
    "New York": "NY",
    "North Carolina": "NC",
    "North Dakota": "ND",
    "Ohio": "OH",
    "Oklahoma": "OK",
    "Oregon": "OR",
    "Pennsylvania": "PA",
    "Rhode Island": "RI",
    "South Carolina": "SC",
    "South Dakota": "SD",
    "Tennessee": "TN",
    "Texas": "TX",
    "Utah": "UT",
    "Vermont": "VT",
    "Virginia": "VA",
    "Washington": "WA",
    "West Virginia": "WV",
    "Wisconsin": "WI",
    "Wyoming": "WY",
    "District of Columbia": "DC",
    "Washington DC": "DC",
    "American Samoa": "AS",
    "Guam": "GU",
    "Northern Mariana Islands": "MP",
    "Puerto Rico": "PR",
    "United States Minor Outlying Islands": "UM",
    "U.S. Virgin Islands": "VI",
    "Virgin Islands, U.S.": "VI",
    "Virgin Islands": "VI",
}

def collapse_metadata( md_loc, output ):
    na_countries = ['Antigua and Barbuda', 'Bahamas', 'Barbados', 'Belize',
                    'Bermuda', 'Canada', 'Costa Rica', 'Cuba',
                    'Dominica', 'Dominican Republic', 'El Salvador', 'Guadeloupe', 'Grenada',
                    'United States', 'Guatemala', 'Haiti', 'Honduras',
                    'Jamaica', 'Mexico', 'Panama', 'Saint Barthélemy', 'Saint Kitts and Nevis',
                    'Saint Lucia', 'Saint Martin', 'Saint Vincent and the Grenadines', 'Sint Maarten']
    md = pd.read_csv( md_loc )
    md["site"] = "Other"
    md.loc[md["country"].isin( na_countries),"site"] = md["country"]
    md.loc[md["country"].isin( ["Canada", "Mexico"] ), "site"] = md["division"]
    md.loc[md["country"] == "United States", "site"] = md["location"]
    md.loc[md["division"].isin( ["Guam", "Puerto Rico", "Virgin Islands", "Northern Mariana Islands"] ), "site"] = md["division"]

    md.loc[(md["country"]=="United States"),"site"] = md.loc[(md["country"]=="United States"),"site"] + "_" + md.loc[(md["country"]=="United States"),"division"].map( us_state_to_abbrev )

    early = md.sort_values( "date_collected" ).head(50)["accession_id"].to_list()
    md.loc[md["accession_id"].isin( early ),"site"] = "root"

    md.to_csv( output, index=False )

if __name__ == "__main__":
    parser = argparse.ArgumentParser( description="Adds column to metadata with most-relevant location data" )

    parser.add_argument( "--input", help="location of metadata", required=True )
    parser.add_argument( "--output", help="location to save appended metadata", required=True )

    args = parser.parse_args()

    collapse_metadata( args.input, args.output )