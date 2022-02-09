rule generate_nulls:
    message: "Shuffle the tips in a tree after grouping by {params.groupby}. Generates null_tree #{wildcards.num}"
    conda: "../envs/general.yaml"
    input:
        tree = config["input_locations"]["tree"],
        metadata = config["input_locations"]["metadata"]
    params:
        groupby = config["generate_nulls"]["groupby"],
        id_col = config["columns"]["id_col"],
        date_col = config["columns"]["date_col"]
    output:
        map_file = temp( "results/temp/null_map_{num}.map" ),
        shuffled_tree = "results/trees/null_{num}.tree"
    script:
        "../scripts/shuffle_tips.py"

rule prune_tree_to_pair:
    message: "Prune {wildcards.status}-{wildcards.num} tree to only sequences from pair: {wildcards.pair}"
    conda: "../envs/general.yaml"
    log: "logs/{pair}.{status}.{num}.pruning.txt"
    group: "pair-split"
    input:
        tree = get_correct_tree,
        metadata = config["input_locations"]["metadata"]
    params:
        pair_list = lambda wildcards: PAIRS[wildcards.pair],
        id_col = config["columns"]["id_col"],
        location_col = config["columns"]["location_col"]
    output:
        pruned_tree = "results/trees/{pair}/{pair}.{status}.{num}.tree"
    script: "../scripts/prune_to_pair.py"

rule compute_phylosor:
    message: "Compute phylosor across time for pair: {wildcards.pair}"
    conda: "../envs/general.yaml"
    log: "logs/{pair}.{status}.{num}.phylosor.txt"
    group: "pair-split"
    input:
        tree = rules.prune_tree_to_pair.output.pruned_tree,
        metadata = config["input_locations"]["metadata"]
    params:
        window_size = config["compute_phylosor"]["window_size"]
    output:
        results = "results/phylosor/{pair}/{pair}.{status}.{num}.csv"
    script: "../scripts/phylosor.R"

rule combine_results:
    message: "Combine phylosor results for all comparisons"
    conda: "../envs/general.yaml"
    log: "logs/combine_results.txt"
    input:
        results = expand( "results/phylosor/{pair}/{pair}.{status}.{num}.csv", pair=PAIRS, status=["actual", "null"], num=range(1,11) )
    output:
        results = "results/output/phylosor_results.csv"
    script: "../scripts/combine_results.py"
