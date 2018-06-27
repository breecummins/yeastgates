import DSGRN

def transform_truth_table(network,truthtables):
    FP_bounds = []
    for tt in truthtables:
        bounds=[]
        for _, t in tt.items():
            bounds.append({network.index(str(k)): t[k] for k in t})
        FP_bounds.append(tuple(bounds))
    return FP_bounds

def truth_table_query(fp_bounds):
    for