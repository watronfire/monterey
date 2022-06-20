def get_correct_tree( wildcards ):
    if wildcards.status == "actual":
        return config["input_locations"]["tree"]
    else:
        return f"results/trees/null_{wildcards.num}.tree"

def get_pair_dict( wildcards ):
    with checkpoints.generate_pairs.get( **wildcards ).output[0].open() as pair_file:
        PAIRS = dict()
        for line in pair_file:
            if not line.startswith( "#" ):
                pair = line.strip().split( "," )
                short_name = line.strip().replace( ",", "-" ).replace( " ", "" )
                PAIRS[short_name] = pair
    return PAIRS

def get_pair_list( wildcards ):
    return get_pair_dict(wildcards)[wildcards.pair]
