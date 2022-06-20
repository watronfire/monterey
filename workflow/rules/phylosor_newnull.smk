rule generate_metadata:
    message: "Remove entries in metadata that aren't found in tree."
    input:
        metadata = config["input_locations"]["metadata"],
        early_chn = config["input_locations"]["chn_metadata"]
    output:
        combined_metadata = "resources/combined_md.csv"
    run:
        import pandas as pd
        ec = pd.read_csv( input.early_chn, parse_dates=["date_collected"] )

        md = pd.read_csv( input.metadata, parse_dates=["date_collected"] )
        md = md.drop( columns="Unnamed: 0" )
        md = pd.concat( [md, ec] )
        md["strain"] = md["strain"].str.slice( start=8 )
        md.to_csv( output.combined_metadata, index=False )

rule metadata_prune:
    message: "Remove tips of tree that aren't found in metadata"
    input:
        tree = config["input_locations"]["tree"],
        metadata = rules.generate_metadata.output.combined_metadata
    output:
        missing = "intermediates/metadata_prune/missing.txt",
        tree = "intermediates/metadata_prune/cog_md.tree",
    run:
        import numpy as np
        md = pd.read_csv( input.metadata, parse_dates=["date_collected"] )

        tree_tips = [i for i in shell( "gotree labels < {input.tree}", iterable=True )]

        # Identify tips that are in tree but aren't in metadata
        missing = np.setdiff1d( tree_tips, md["strain"].to_list() )
        print( f"{len(missing)} of {len(md['strain'])} tips missing" )

        with open( output.missing, "w" ) as missing_file:
            missing_file.write( "\n".join( missing ) )

        # Prune tips that are missing
        shell( "gotree prune --tipfile {output.missing} < {input.tree} > {output.tree}" )


rule rename_tree_to_accession:
    message: "Rename tips of tree to accession IDs from strain names"
    input:
        tree = rules.metadata_prune.output.tree,
        metadata = rules.generate_metadata.output.combined_metadata
    output:
        renames = "intermediates/rename_tree/renames.tsv",
        tree = "intermediates/rename_tree/cog_accession.tree"
    shell:
        """
        cut -f1,2 -d, combined_md.csv | sed "s/,/\t/g" > {output.renames} &&
        gotree rename --map {output.renames} < {input.tree} > {output.tree}
        """

rule generate_pairs:
    message: "Do a robust search through the metadata for all locations which have greater than {params.sequences} and greater than {params.completeness} epiweeks covered."
    conda: "../envs/general.yaml"
    log: "logs/generate_pairs.txt"
    input:
        metadata = rules.generate_metadata.output.combined_metadata
    params:
        sequences = config["pairs"]["min_sequences"],
        completeness = config["pairs"]["min_completeness"]
    output:
        pairs = "intermediates/pairs/pairs.txt",
        pair_graph = "results/reports/pair_graph.pdf"
    shell:
        """
        python workflow/scripts/generate_pairs.py \
            --metadata {input.metadata} \
            --min-sequences {params.sequences} \
            --min-completeness {params.completeness} \
            --output {output.pairs} \
            --graph {output.pair_graph}
        """

rule prune_tree_to_pair_newnull:
    message: "Prune tree to only sequences from pair: {wildcards.pair}"
    conda: "../envs/general.yaml"
    log: "logs/{pair}.pruning.txt"
    input:
        tree = rules.metadata_prune.output.tree,
        metadata = rules.generate_metadata.output.combined_metadata
    params:
        pair_list = lambda wildcards: PAIRS[wildcards.pair],
        id_col = config["columns"]["id_col"],
        date_col = config["columns"]["date_col"],
        location_col = config["columns"]["location_col"],
        date_max = config["prune_tree_to_pair"]["date_max"]
    output:
        pruned_tree = "results/trees/{pair}/{pair}.tree"
    shell:
        """
        python workflow/scripts/prune_to_pair.py \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --id {params.id_col} \
            --date-col {params.date_col} \
            --location-col {params.location_col} \
            --pair {params.pair_list:q} \
            --max-date {params.date_max} \
            --output {output.pruned_tree}
        """

rule compute_phylosor_newnull:
    message: "Compute {wildcards.status} phylosor across time for pair: {wildcards.pair}"
    conda: "../envs/general.yaml"
    log: "logs/{pair}.{status}.{num}.phylosor_newnull.txt"
    input:
        tree = rules.prune_tree_to_pair_newnull.output.pruned_tree,
        metadata = rules.generate_metadata.output.combined_metadata
    params:
        pair_list = lambda wildcards: PAIRS[wildcards.pair],
        window_size = config["compute_phylosor"]["window_size"],
        shuffle = lambda wildcards: "--shuffle" if wildcards.status == "null" else ""
    output:
        results = "results/phylosor_newnull/{pair}/{pair}.{status}.{num}.csv"
    shell:
        """
        python workflow/scripts/phylosor_table.py \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --pair-list {params.pair_list:q} \
            --window-size {params.window_size} \
            {params.shuffle} \
            --output {output.results} \
        """

rule compute_hill:
    message: "Compute {wildcards.status} hill number across time for pair: {wildcards.pair}"
    conda: "../envs/general.yaml"
    log: "logs/{pair}.{status}.{num}.hill.txt"
    input:
        tree = rules.prune_tree_to_pair_newnull.output.pruned_tree,
        metadata = config["input_locations"]["metadata"]
    params:
        pair1 = lambda wildcards: PAIRS[wildcards.pair][0],
        pair2 = lambda wildcards: PAIRS[wildcards.pair][1],
        window_size = config["compute_phylosor"]["window_size"],
        shuffle = lambda wildcards: "--shuffle" if wildcards.status == "null" else ""
    output:
        results = "results/hill/{pair}/{pair}.{status}.{num}.csv"
    shell:
        """
        Rscript workflow/scripts/hill_monthly.R \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --pair1 {params.pair1:q} \
            --pair2 {params.pair2:q} \
            --window-size {params.window_size} \
            {params.shuffle} \
            --output {output.results} \
        """

#rule combine_hill:
#    message: "Combine hill results for all comparisons"
#    conda: "../envs/general.yaml"
#    log: "logs/combine_results.txt"
#    input:
#        results_nulls = expand( "results/hill/{pair}/{pair}.null.{num}.csv",pair=JS_PAIRS,num=range( 1,11 ) ),
#        results_actual = expand( "results/hill/{pair}/{pair}.actual.{num}.csv",pair=JS_PAIRS,num=[1] )
#    output:
#        results = "results/output/hill_results.csv"
#    shell:
#        """
#        python workflow/scripts/combine_results.py \
#            {input.results_actual} \
#            {input.results_nulls} \
#            {output.results}
#        """
#
#rule combine_results_newnull:
#    message: "Combine phylosor results for all comparisons"
#    conda: "../envs/general.yaml"
#    log: "logs/combine_results.txt"
#    input:
#        results_nulls = expand( "results/phylosor_newnull/{pair}/{pair}.null.{num}.csv",pair=PAIRS,num=range( 1,11 ) ),
#        results_actual = expand( "results/phylosor_newnull/{pair}/{pair}.actual.{num}.csv",pair=PAIRS,num=[1] )
#    output:
#        results = "results/output/phylosor_newnull_results.csv"
#    shell:
#        """
#        python workflow/scripts/combine_results.py \
#            {input.results_actual} \
#            {input.results_nulls} \
#            {output.results}
#        """
#
#rule plot_results_newnull:
#    message: "Plot phylosor metric for pair: {wildcards.pair}"
#    conda: "../envs/general.yaml"
#    log: "logs/{pair}.plotting.log"
#    input:
#        metadata = config["input_locations"]["metadata"],
#        results = rules.combine_results_newnull.output.results
#    params:
#        pair_list = lambda wildcards: PAIRS[wildcards.pair]
#    output:
#        plot = "results/phylosor_plot/{pair}.phylosor.pdf"
#    script:
#        "../scripts/plot_phylosor.py"
