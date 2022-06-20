def get_correct_tree( wildcards ):
    if wildcards.status == "actual":
        return config["input_locations"]["tree"]
    else:
        return f"results/trees/null_{wildcards.num}.tree"

def get_pair_dict( wildcards, cp ):
    with cp.generate_pairs.get( **wildcards ).output[0].open() as pair_file:
        PAIRS = dict()
        for line in pair_file:
            if not line.startswith( "#" ):
                pair = line.strip().split( "," )
                short_name = line.strip().replace( ",", "-" ).replace( " ", "" )
                PAIRS[short_name] = pair
    return PAIRS

#def get_pair_list( wildcards ):
#    return get_pair_dict(wildcards)[wildcards.pair]

def generate_phylosor_results( wildcards ):
    global checkpoints
    PAIRS = get_pair_dict( wildcards, checkpoints )

    results = expand( "results/phylosor_newnull/{pair}/{pair}.null.{num}.csv", pair=PAIRS, num=range( 1,11 ) )
    results.extend( expand( "results/phylosor_newnull/{pair}/{pair}.actual.{num}.csv", pair=get_pair_dict, num=[1] ) )
    return results