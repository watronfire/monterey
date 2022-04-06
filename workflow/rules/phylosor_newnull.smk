rule prune_tree_to_pair_newnull:
    message: "Prune tree to only sequences from pair: {wildcards.pair}"
    conda: "../envs/general.yaml"
    log: "logs/{pair}.pruning.txt"
    input:
        tree = config["input_locations"]["tree"],
        metadata = config["input_locations"]["metadata"]
    params:
        pair_list = lambda wildcards: PAIRS[wildcards.pair],
        id_col = config["columns"]["id_col"],
        date_col = config["columns"]["date_col"],
        location_col = config["columns"]["location_col"],
        date_max = config["prune_tree_to_pair"]["date_max"]
    output:
        pruned_tree = "results/trees/{pair}/{pair}.tree"
    script: "../scripts/prune_to_pair.py"

rule compute_phylosor_newnull:
    message: "Compute {wildcards.status} phylosor across time for pair: {wildcards.pair}"
    conda: "../envs/general.yaml"
    log: "logs/{pair}.{status}.{num}.phylosor_newnull.txt"
    input:
        tree = rules.prune_tree_to_pair_newnull.output.pruned_tree,
        metadata = config["input_locations"]["metadata"]
    params:
        pair_list = lambda wildcards: PAIRS[wildcards.pair],
        window_size = config["compute_phylosor"]["window_size"],
        shuffle = lambda wildcards: wildcards.status == "null"
    output:
        results = "results/phylosor_newnull/{pair}/{pair}.{status}.{num}.csv"
    script: "../scripts/phylosor_table.py"

rule combine_results_newnull:
    message: "Combine phylosor results for all comparisons"
    conda: "../envs/general.yaml"
    log: "logs/combine_results.txt"
    input:
        results_nulls = expand( "results/phylosor_newnull/{pair}/{pair}.null.{num}.csv",pair=PAIRS,num=range( 1,11 ) ),
        results_actual = expand( "results/phylosor_newnull/{pair}/{pair}.actual.{num}.csv",pair=PAIRS,num=[1] )
    output:
        results = "results/output/phylosor_newnull_results.csv"
    script: "../scripts/combine_results.py"

rule plot_results_newnull:
    message: "Plot phylosor metric for pair: {wildcards.pair}"
    conda: "../env/general.yaml"
    log: "logs/{pair}.plotting.log"
    input:
        metadata = config["input_locations"]["metadata"],
        results = rules.combine_results_newnull.output.results
    params:
        pair_list = lambda wildcards: PAIRS[wildcards.pair]
    output:
        plot = "results/phylosor_plot/{pair}.phylosor.pdf"
    script:
        "../scripts/plot_phylosor.py"