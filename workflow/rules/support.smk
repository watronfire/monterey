
def get_correct_tree( wildcards ):
    if wildcards.status == "actual":
        return config["input_locations"]["tree"]
    else:
        return f"results/trees/null_{wildcards.num}.tree"