import pandas as pd
import time
from epiweeks import Week
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

def timeseries_formatting( axis, spines=["bottom"], which="y", title=None, title_loc="left", ylabel=None, xlabel=None, xsize=12, ysize=12, xlims=None, ylims=None ):
    # Properly label timeseries
    axis.xaxis.set_minor_locator( mdates.MonthLocator() )
    axis.xaxis.set_minor_formatter( mdates.DateFormatter( '%b' ) )
    axis.xaxis.set_major_locator( mdates.YearLocator() )
    axis.xaxis.set_major_formatter( mdates.DateFormatter( '%Y %b' ) )

    # Remove spines
    [axis.spines[j].set_visible( False ) for j in axis.spines if j not in spines]

    # Format axes ticks
    axis.tick_params( axis="x", bottom=False, which="both", rotation=90, labelbottom=True, labelsize=xsize )
    axis.tick_params( axis="y", left=False, which="both", labelleft=True, labelsize=ysize )

    # Label axes
    axis.set_xlabel( xlabel, fontsize=xsize )
    axis.set_ylabel( ylabel, fontsize=ysize )
    axis.set_title( title, loc=title_loc )

    # Add a simple grid
    axis.grid( which="both", axis=which, linewidth=1, color="#F1F1F1", zorder=1 )

    # Add the xlims
    if xlims:
        axis.set_xlim( xlims )
    if ylims:
        axis.set_ylim( ylims )


def plot_sampling( axis, df, colors=["#56B4E9", "#118f60"], order=None, fixed=False ):
    plot_df = df.pivot_table( index="week", columns="site", values="pangolin_lineage", aggfunc="count" ).fillna(0.0)

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

    timeseries_formatting( axis, which="both", )
    timeseries_formatting( axis2, which="both", spines=["top"] )
    axis.tick_params( which="both", labelbottom=False )
    axis.set_ylabel( plots[0].name )
    axis2.set_ylabel( plots[1].name )


def plot_jensenshannon( axis, res ):
    axis.plot( res.index, res["js_mean"], color="black", zorder=10 )
    axis.fill_between( x=res.index, y1=res["js_5"], y2=res["js_95"], color="black", alpha=0.1, linewidth=0, zorder=9 )
    timeseries_formatting( axis, ylabel="Jensen-Shannon Distance", which="both", ylims=(0,1) )


def load_metadata( md_loc, pair, verbose=True ):
    if verbose:
        print( "Loading metadata...", end="" )
    starting_time = time.time()
    metadata = pd.read_csv( md_loc, usecols=["accession_id", "date_collected", "pangolin_lineage", "site"], parse_dates=["date_collected"] )
    metadata = metadata.loc[metadata["site"].isin(pair)]
    metadata["week"] = metadata["date_collected"].apply( lambda x: Week.fromdate( x ).startdate() )
    metadata["month"] = metadata["date_collected"].to_numpy().astype('datetime64[M]')
    if verbose:
        print( f"Found {metadata.shape[0]} genomes in {time.time() - starting_time:.1f} seconds" )
    return metadata


def genome_bootstrap( df, fun, replicates=100 ):
    result_list = list()
    for _ in range( replicates ):
        trial = df.groupby( "site" ).sample( frac=1, replace=True )
        result_list.append( fun( trial ) )
    res = np.quantile( result_list, q=[0.025, 0.5, 0.975] )
    return {k:j for k,j in zip( ["js_5", "js_mean", "js_95"], res )}


def apply_jensenshannon( df ):
    evaluate = df.pivot_table( index="pangolin_lineage", columns="site", values="accession_id", aggfunc="count", fill_value=0 )
    evaluate = evaluate.apply( lambda x: x / x.sum() )
    try:
        return jensenshannon( evaluate.iloc[:,0], evaluate.iloc[:,1], base=2 )
    except IndexError:
        return np.nan


if __name__ == "__main__":
    md = load_metadata( snakemake.input.metadata, snakemake.params.pair_list )
    results = md.groupby( snakemake.params.resolution ).apply( genome_bootstrap, fun=apply_jensenshannon )
    results = pd.DataFrame( results.to_list(), index=results.index )
    results = results.dropna()
    results["siteA"] = snakemake.params.pair_list[0]
    results["siteB"] = snakemake.params.pair_list[1]
    results.to_csv( snakemake.output.results )

    # Plot stuff
    fig, ax = plt.subplots( dpi=200, figsize=(10,8), nrows=2, sharex=True )
    plot_sampling( ax[0], md, order=snakemake.params.pair_list )
    plot_jensenshannon( ax[1], results )
    ax[0].set_title( "Sequence counts", loc="left" )
    ax[1].set_title( "PANGO lineage dissimilarity", loc="left" )
    plt.tight_layout()
    plt.savefig( snakemake.output.plot )