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
        tree = rules.prune_tree_to_pair.output.pruned_tree,
        metadata = config["input_locations"]["metadata"]
    params:
        pair_list = lambda wildcards: PAIRS[wildcards.pair],
        window_size = config["compute_phylosor"]["window_size"],
        shuffle = lambda wildcards: wildcards.status == "null"
    output:
        results = "results/phylosor_newnull/{pair}/{pair}.{status}.{num}.csv"
    script: "../scripts/phylosor_table.py"

