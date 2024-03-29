import argparse
import pandas as pd
from epiweeks import Week

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

def get_passing_sites( df, min_sequences=1000, min_completeness=0.75 ):
    # sequences per site
    seqs = df["site"].value_counts()

    # fraction of sequences per week per site.
    completeness = df.pivot_table( columns="site", index="week", values="date_collected", aggfunc="count" )
    completeness = completeness.loc[completeness.index > pd.to_datetime( "2020-02-23" )]
    completeness = completeness.count() / completeness.shape[0]

    requirements = pd.concat( [seqs, completeness], axis=1, ignore_index=False )
    requirements.columns = ["sequences", "completeness"]

    return requirements.loc[(requirements["sequences"]>=min_sequences)&(requirements["completeness"]>=min_completeness)].index.to_list()

def collapse_metadata( md_loc, output, min_sequences, min_completeness ):
    md = pd.read_csv( md_loc, parse_dates=["date_collected"] )
    md["week"] = md["date_collected"].apply( lambda x: Week.fromdate( x ).startdate() )

    na = md.loc[md["country"].isin( ["United States", "Mexico", "Canada"] )].copy()
    na.loc[na["country"] == "Canada", "site"] = na["division"] + "_CAN"
    na.loc[na["country"] == "Mexico", "site"] = na["division"] + "_MEX"
    us = na.loc[(na["country"] == "United States") & ~na["division"].isna()].copy()
    us["site"] = us["location"] + "_" + us["division"].map( us_state_to_abbrev )

    passing = get_passing_sites( us, min_sequences=min_sequences, min_completeness=min_completeness )
    us.loc[~us["site"].isin( passing ), "site"] = us["division"] + "_USA"

    md_add = pd.concat(
        [md.loc[~md["country"].isin( ["United States", "Mexico", "Canada"] )], na.loc[na["country"] != "United States"],
         us], ignore_index=True )
    md_add["site"] = md_add["site"].fillna( "Other" )
    early = md_add.sort_values( "date_collected" ).head( 50 )["accession_id"].to_list()
    md_add.loc[md_add["accession_id"].isin( early ), "site"] = "root"

    md_add.to_csv( output, index=False )


def collapse_metadata_alternative( md_loc, output ):
    na_countries = ['Antigua and Barbuda', 'Bahamas', 'Barbados', 'Belize',
                    'Bermuda', 'Canada', 'Costa Rica', 'Cuba',
                    'Dominica', 'Dominican Republic', 'El Salvador', 'Guadeloupe', 'Grenada',
                    'United States', 'Guatemala', 'Haiti', 'Honduras',
                    'Jamaica', 'Mexico', 'Panama', 'Saint Barthélemy', 'Saint Kitts and Nevis',
                    'Saint Lucia', 'Saint Martin', 'Saint Vincent and the Grenadines', 'Sint Maarten']
    md = pd.read_csv( md_loc )
    md["site"] = "Other"
    md.loc[md["country"].isin( na_countries),"site"] = md["country"]
    md.loc[md["country"].isin( ["Canada", "Mexico", "United States"] ), "site"] = md["division"]
    md.loc[md["division"] == "California", "site"] = md["location"]
    #md.loc[md["division"].isin( ["Guam", "Puerto Rico", "Virgin Islands", "Northern Mariana Islands"] ), "site"] = md["division"]

    md.loc[(md["country"]=="United States")&(md["division"]!="California"),"site"] = md.loc[(md["country"]=="United States")&(md["division"]!="California"),"site"] + "_USA"
    md.loc[(md["country"]=="Canada"), "site"] = md.loc[(md["country"]=="Canada"), "site"] + "_CAN"
    md.loc[(md["country"]=="Mexico"), "site"] = md.loc[(md["country"]=="Mexico"), "site"] + "_MEX"
    md.loc[md["division"]=="California", "site"] = md.loc[md["division"]=="California", "site"] + "_CA"

    early = md.sort_values( "date_collected" ).head(50)["accession_id"].to_list()
    md.loc[md["accession_id"].isin( early ),"site"] = "root"

    md.to_csv( output, index=False )


if __name__ == "__main__":
    parser = argparse.ArgumentParser( description="Adds column to metadata with most-relevant location data" )

    parser.add_argument( "--input", help="location of metadata", required=True )
    parser.add_argument( "--output", help="location to save appended metadata", required=True )
    parser.add_argument( "--min-sequences", help="Keep locations with at least this many sequences", type=int, required=True )
    parser.add_argument( "--min-completeness", help="Keep locations with sequences collected from at least this many epiweeks", type=float, required=True )
    parser.add_argument( "--alternative", action="store_true", help="Use non-greedy collapsing algorithm" )

    args = parser.parse_args()

    if args.alternative:
        collapse_metadata_alternative( args.input, args.output )
    else:
        collapse_metadata( args.input, args.output, args.min_sequences, args.min_completeness )