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

    return requirements.loc[(requirements["sequences"]>=min_completeness)&(requirements["completeness"]>=min_completeness)].index.to_list()

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

if __name__ == "__main__":
    parser = argparse.ArgumentParser( description="Adds column to metadata with most-relevant location data" )

    parser.add_argument( "--input", help="location of metadata", required=True )
    parser.add_argument( "--output", help="location to save appended metadata", required=True )
    parser.add_argument( "--min-sequences", help="Keep locations with at least this many sequences", required=True )
    parser.add_argument( "--min-completeness", help="Keep locations with sequences collected from at least this many epiweeks", required=True )

    args = parser.parse_args()

    collapse_metadata( args.input, args.output )