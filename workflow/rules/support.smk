def get_correct_tree( wildcards ):
    if wildcards.status == "actual":
        return config["input_locations"]["tree"]
    else:
        return f"results/trees/null_{wildcards.num}.tree"

def get_pair_dict( wildcards ):
    global checkpoints
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

def generate_phylosor_results( wildcards ):
    PAIRS = get_pair_dict( wildcards )

    #results = expand( "results/phylosor/{pair}/{pair}.null.{num}.csv", pair=PAIRS, num=range( 1,11 ) )
    results = expand( "results/phylosor/{kind}.{num}/{pair}/{pair}.actual.csv", pair=PAIRS, kind=["actual"], num=[1] )
    results.extende( expand( "results/phylosor/{kind}.{num}/{pair}/{pair}.actual.csv", pair=PAIRS, kind=["fraction", "count"], num=range(1, 11) ) )
    return results

def generate_hill_results( wildcards ):
    PAIRS = get_pair_dict( wildcards )

    results = expand( "results/hill/{pair}/{pair}.null.{num}.csv", pair=PAIRS, num=range( 1,11 ) )
    results.extend( expand( "results/hill/{pair}/{pair}.actual.{num}.csv", pair=PAIRS, num=[1] ) )
    return results
