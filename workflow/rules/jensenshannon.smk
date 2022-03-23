rule compute_jensenshannon:
    message: "Compute Jensen-Shannon distance over time for pair: {wildcards.pair}"
    conda: "../envs/general.yaml"
    log: "logs/{pair}.js.txt"
    input:
        metadata = config["input_locations"]["metadata"]
    params:
        pair_list = lambda wildcards: JS_PAIRS[wildcards.pair],
        resolution = config["compute_js"]["resolution"]
    output:
        results = "results/js/{pair}.js.csv",
        plot = "results/js/{pair}.js.pdf"
    script: "../scripts/jensen_shannon_distance.py"
