import pandas as pd
import time
from epiweeks import Week
from matplotlib.lines import Line2D
from scipy.spatial.distance import jensenshannon
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.dates as mdates
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl

prop = mpl.font_manager.FontProperties('Roboto')
mpl.rcParams['font.sans-serif'] = prop.get_name()
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.weight']=300
mpl.rcParams['axes.labelweight']=300
mpl.rcParams['font.size']=16

mpl.rcParams['figure.dpi'] = 200

COLOR = '#343434'
mpl.rcParams['text.color'] = COLOR
mpl.rcParams['axes.labelcolor'] = COLOR
mpl.rcParams['xtick.color'] = COLOR
mpl.rcParams['ytick.color'] = COLOR

def timeseries_formatting( ax, spines=["bottom"], which="y", title=None, ylabel=None, xlabel=None, xsize=10, ysize=10, xlims=None, ylims=None ):
    # Properly label timeseries
    ax.xaxis.set_minor_locator( mdates.MonthLocator() )
    ax.xaxis.set_minor_formatter( mdates.DateFormatter( '%b' ) )
    ax.xaxis.set_major_locator( mdates.YearLocator() )
    ax.xaxis.set_major_formatter( mdates.DateFormatter( '%Y %b' ) )

    # Remove spines
    [ax.spines[j].set_visible(False) for j in ax.spines if j not in spines]

    # Format axes ticks
    ax.tick_params( axis="x", bottom=False, which="both", rotation=90, labelbottom=True, labelsize=xsize )
    ax.tick_params( axis="y", left=False, which="both", labelleft=True, labelsize=ysize )

    # Label axes
    ax.set_xlabel( xlabel, fontsize=xsize )
    ax.set_ylabel( ylabel, fontsize=ysize )
    ax.set_title( title)

    # Add a simple grid
    ax.grid( which="both", axis=which, linewidth=1, color="#F1F1F1", zorder=1 )

    # Add the xlims
    if xlims:
        ax.set_xlim( xlims )
    if ylims:
        ax.set_ylim( ylims )

def plot_sampling( axis, df, colors=["#E69F00", "#56B4E9"], order=None, fixed=False ):
    plot_df = df.pivot_table( index="epiweek", columns="site", values="strain", aggfunc="count" ).fillna(0.0)

    if order:
        plots = [plot_df[i] for i in order]
    else:
        plots = [plot_df.iloc[:,0], plot_df.iloc[:,1]]

    axis.bar( x=plot_df.index, height=plots[0], color=colors[0], width=6, zorder=10 )

    divider = make_axes_locatable( axis )
    axis2 = divider.append_axes( "bottom", size="100%", pad=0, sharex=axis )
    axis2.bar( x=plot_df.index, height=plots[1], color=colors[1], width=6, zorder=10 )

    if fixed:
        if axis.get_ylim() < axis2.get_ylim():
            axis.set_ylim( axis2.get_ylim() )
        else:
            axis2.set_ylim( axis.get_ylim() )

    axis2.invert_yaxis()

    timeseries_formatting( axis )
    timeseries_formatting( axis2, spines=["top"] )
    axis.tick_params( which="both", labelbottom=False )
    axis.set_ylabel( plots[0].name )
    axis2.set_ylabel( plots[1].name )

def plot_within_sampling( axis, df, name, color="#56B4E9",  ):
    plot_df = df["epiweek"].value_counts().reset_index().sort_index()
    axis.bar( x="index", height="epiweek", data=plot_df, color=color, width=6, zorder=10 )
    timeseries_formatting( axis, spines=["bottom"] )
    axis.set_ylabel( name )

def plot_phylosor_nulls( axis, df, focus, color, col="value", missing=True, title=False ):
    plot_df = format_for_focus( df, focus[0] )
    plot_df = plot_df.loc[plot_df["from"]==focus[1]].dropna()

    for i in sorted( plot_df["kind"].unique(), reverse=True ):
        axis.plot( "date", col, color="#DBDBDB" if i.startswith( "null" ) else color, data=plot_df.loc[plot_df["kind"]==i] )

    legend = [Line2D([0], [0], linestyle='none', marker='o', color=color, label="True locations", markersize=12 ),
              Line2D([0], [0], linestyle='none', marker='o', color="#DBDBDB", label="Shuffled locations", markersize=12 )]
    axis.legend( loc="upper left", handletextpad=0.1, handles=legend, frameon=False, fontsize=12 )

    if title:
        axis.set_title( " and ".join( focus ), loc="left" )

def plot_phylosor( axis, df, focus, color, missing=True, normalized=False ):

    plotting_values = { False: "value",
                        "z" : "corrected_z",
                        "sub" : "corrected_sub",
                        "div" : "corrected_div" }
    pvalue = plotting_values[normalized]

    plot_df = format_for_focus( df, focus[0] )
    plot_df = plot_df.loc[(plot_df["to"]==focus[0])&(plot_df["from"]==focus[1])].dropna()
    if not normalized:
        plot_df = plot_df.loc[plot_df["kind"]=="actual1"]

    if not missing:
        plot_df.loc[plot_df["value"]==0,"value"] = np.nan

    axis.plot( "date", pvalue, color=color, data=plot_df, zorder=10 )
    axis.fill_between( x="date", y1="corrected_upper", y2="corrected_lower", data=plot_df, color="black", alpha=0.14, linewidth=0, zorder=9 )

def format_for_focus( input_df, focus ):
    df = input_df.copy()
    df["to"] = df["siteA"]
    df["from"] = df["siteB"]
    df.loc[df["siteB"]==focus,"to"] = df["siteB"]
    df.loc[df["siteB"]==focus,"from"] = df["siteA"]
    df = df.loc[df["to"]==focus]
    return df

def load_metadata( md_loc, pair ) :
    return_df = pd.read_csv( md_loc, usecols=["date_collected", "site", "strain"], parse_dates=["date_collected"] )
    return_df = return_df.loc[return_df["site"].isin( pair )]
    return_df["epiweek"] = return_df["date_collected"].apply( lambda x: Week.fromdate( x ).startdate() )
    return return_df

def get_corrected_df( entry, col="value" ):
    a = entry.loc[entry["kind"].str.startswith("actual")]
    n = entry.loc[entry["kind"].str.startswith("null")]
    n = n.groupby( ["siteA", "siteB", "date"] )[col].agg(
        null_upper=lambda x: x.quantile( 0.975 ),
        null_lower=lambda x: x.quantile( 0.025 ),
        null_mean="mean",
        null_std="std",
        null_count="count" )
    a = a.merge( n, left_on=["siteA", "siteB", "date"], right_index=True, validate="1:1" )
    a["corrected_sub"] = a["null_mean"] - a[col]
    a["corrected_upper"] = a["null_lower"] - a[col]
    a["corrected_lower"] = a["null_upper"] - a[col]
    return a

def load_results( results_loc, pair ) :
    return_df = pd.read_csv( results_loc, parse_dates=["date"] )
    return_df["kind"] = return_df["kind"].fillna( "null" )
    return_df["kind"] = return_df["kind"] + return_df["num"].astype( str )
    return_df = return_df.drop( columns=["num"] )
    return_df = return_df.loc[(return_df["siteA"]==pair[0])&(return_df["siteB"]==pair[1])]

    df_corrected = get_corrected_df( return_df )

    return return_df, df_corrected

def plot_phylosor_all( md_loc, results_loc, pair_list, output ):
    md = load_metadata( md_loc, pair_list )
    results, results_corr = load_results( results_loc, pair_list )

    fig, ax = plt.subplots( dpi=200, figsize=(10,13), nrows=4, sharex=True )
    if pair_list[0] == pair_list[1]:
        plot_within_sampling( ax[0], md, pair_list[0] )
    else:    
        plot_sampling( ax[0], md, order=pair_list )
    plot_phylosor_nulls( ax[1], results, pair_list, COLOR )
    timeseries_formatting( ax[1], ylabel="Proportion of\nbranch length shared", ylims=[0,1] )
    plot_phylosor( ax[2], results_corr, pair_list, COLOR, normalized="sub" )
    timeseries_formatting( ax[2], ylabel="Difference to mixed-model", ylims=[0,1] )
    plot_phylosor_nulls( ax[3], results, pair_list, COLOR, col="value_turn" )
    timeseries_formatting( ax[3], ylabel="$PhyloSor_{Turn}$", ylims=[0,1] )
    plt.suptitle( " <-> ".join( pair_list ), fontsize=12 )
    plt.tight_layout()
    plt.savefig( output )

if __name__ == "__main__":
    plot_phylosor_all(
        snakemake.input.metadata,
        snakemake.input.results,
        snakemake.params.pair_list,
        snakemake.output.plot,
    )
