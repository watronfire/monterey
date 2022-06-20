def get_correct_tree( wildcards ):
    if wildcards.status == "actual":
        return config["input_locations"]["tree"]
    else:
        return f"results/trees/null_{wildcards.num}.tree"

def get_pair_dict():
    with checkpoints.generate_pairs.get().output[0].open() as pair_file:
        PAIRS = dict()
        for line in pair_file:
            if not line.startswith( "#" ):
                pair = line.strip().split( "," )
                short_name = line.strip().replace( ",", "-" ).replace( " ", "" )
                PAIRS[short_name] = pair

def get_pair_list( wildcards ):
    return get_pair_dict()[wildcards.pair]