rule generate_nulls:
    message: "Shuffle the tips in a tree after grouping by {params.groupby}. Generates null_tree #{wildcards.null}"
    conda: "../envs/general.yaml"
    input:
        tree = config["input_locations"]["tree"],
        metadata = config["input_locations"]["metadata"]
    params:
        groupby = config["generate_nulls"]["groupby"],
        id_col = config["generate_nulls"]["id_col"],
        date_col = config["generate_nulls"]["date_col"]
    output:
        map_file = temp( "results/temp/null_map_{null}.map" ),
        shuffled_tree = "results/trees/null_{null}.tree"
    script:
        "../scripts/shuffle_tips.py"