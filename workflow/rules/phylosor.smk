rule generate_metadata:
    message: "Remove entries in metadata that aren't found in tree."
    input:
        metadata = config["input_locations"]["metadata"],
        early_chn = config["input_locations"]["chn_metadata"]
    output:
        combined_metadata = "resources/combined_md.csv.gz"
    run:
        import pandas as pd
        ec = pd.read_csv( input.early_chn, parse_dates=["date_collected"] )

        md = pd.read_csv( input.metadata, parse_dates=["date_collected"] )
        if "Unnamed: 0" in md.columns:
            md = md.drop( columns="Unnamed: 0" )
        md = pd.concat( [md, ec] )
        md["strain"] = md["strain"].str.slice( start=8 )
        md.to_csv( output.combined_metadata, index=False )

rule metadata_prune:
    message: "Remove tips of tree that aren't found in metadata"
    input:
        tree = config["input_locations"]["tree"],
        metadata = rules.generate_metadata.output.combined_metadata,
        sampling = config["metadata_prune"]["config"]
    params:
        kind = lambda wildcards: wildcards.kind
    output:
        missing = "intermediates/metadata_prune/missing.{kind}.{num}.txt",
        tree = "intermediates/metadata_prune/cog_md.{kind}.{num}.tree",
    shell:
        """
        python metadata_prune \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --sampling {input.sampling} \
            --location "San Diego" \
            --kind {params.kind} \
            --missing {output.missing} \
            --output {output.tree} \
        """


#rule rename_tree_to_accession:
#    message: "Rename tips of tree to accession IDs from strain names"
#    input:
#        tree = rules.metadata_prune.output.tree,
#        metadata = rules.generate_metadata.output.combined_metadata
#    output:
#        renames = "intermediates/rename_tree/renames.tsv",
#        tree = "intermediates/rename_tree/cog_accession.tree"
#    shell:
#        """
#        cut -f1,2 -d, {input.metadata} | sed "s/,/\t/g" > {output.renames} &&
#        gotree rename --map {output.renames} < {input.tree} > {output.tree}
#        """

rule collapse_location_in_metadata:
    message: "Add column to metadata with most-relevant location data"
    input:
        metadata = rules.generate_metadata.output.combined_metadata
    params:
        sequences = config["pairs"]["min_sequences"],
        completeness = config["pairs"]["min_completeness"],
        alternative = "--alternative" if config["collapse_locations"]["sd_focus"] else ""
    output:
        collapsed_metadata = "intermediates/collapse_location/metadata.csv"
    shell:
        """
        python workflow/scripts/collapse_location.py \
            --input {input.metadata} \
            --min-sequences {params.sequences} \
            --min-completeness {params.completeness} \
            {params.alternative} \
            --output {output.collapsed_metadata}
        """

checkpoint generate_pairs:
    message: "Do a robust search through the metadata for all locations which have greater than {params.sequences} and greater than {params.completeness} epiweeks covered."
    conda: "../envs/general.yaml"
    log: "logs/generate_pairs.txt"
    input:
        metadata = rules.collapse_location_in_metadata.output.collapsed_metadata
    params:
        sequences = config["pairs"]["min_sequences"],
        completeness = config["pairs"]["min_completeness"],
        focus = "--locations 'San Diego_CA'" if config["collapse_locations"]["sd_focus"] else ""
    output:
        pairs = "intermediates/pairs/pairs.txt",
        pair_graph = "results/reports/pair_graph.pdf",
        summary = "intermediates/pairs/pairs.csv"
    shell:
        """
        python workflow/scripts/generate_pairs.py \
            --metadata {input.metadata} \
            --min-sequences {params.sequences} \
            --min-completeness {params.completeness} \
            --output {output.pairs} \
            --graph {output.pair_graph} \
            {params.focus} \
            --summary {output.summary}
       """

rule prune_tree_to_pair:
    message: "Prune tree to only sequences from pair: {wildcards.pair}"
    conda: "../envs/general.yaml"
    input:
        tree = rules.metadata_prune.output.tree,
        metadata = rules.collapse_location_in_metadata.output.collapsed_metadata
    params:
        pair_list = get_pair_list,
        id_col = config["columns"]["id_col"],
        date_col = config["columns"]["date_col"],
        location_col = config["columns"]["location_col"],
        date_max = config["prune_tree_to_pair"]["date_max"]
    output:
        pruned_tree = "results/trees/{kind}.{num}/{pair}/{pair}.tree"
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

rule compute_phylosor:
    message: "Compute {wildcards.status} phylosor across time for pair: {wildcards.pair}"
    conda: "../envs/general.yaml"
    input:
        tree = rules.prune_tree_to_pair.output.pruned_tree,
        metadata = rules.collapse_location_in_metadata.output.collapsed_metadata
    params:
        pair_list = get_pair_list,
        window_size = config["compute_phylosor"]["window_size"],
        shuffle = lambda wildcards: "--shuffle" if wildcards.status == "null" else ""
    output:
        results = "results/phylosor/{kind}.{num}/{pair}/{pair}.{status}.csv"
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


rule combine_results:
    message: "Combine phylosor results for all comparisons"
    conda: "../envs/general.yaml"
    log: "logs/combine_results.txt"
    input:
        results = generate_phylosor_results
    output:
        results = "results/output/phylosor_results.csv"
    shell:
        """
        python workflow/scripts/combine_results.py \
            {input.results} \
            {output.results}
        """
