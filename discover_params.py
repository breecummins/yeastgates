import DSGRN,json
from truth_table_fixed_pts import do_all_queries

in1 = "s1"
in2 = "s2"
out = "GFP"
topvalout = 1

# # NAND
# dbname = "NAND.db"
# network = "NAND.txt"
# topval1 = 2
# topval2 = 2
# print_output=True
# design_spec = ["C","C","2","2","8","2"]

# OR
dbname = "OR.db"
network = "OR.txt"
topval1 = 2
topval2 = 2
print_output=True
design_spec = ["C","C","8","2"]

def get_one_change_params():
    sfile, pfile = do_all_queries(dbname, in1, topval1, in2, topval2, out, topvalout, print_output)

    summary = json.load(open(sfile))
    params = json.load(open(pfile))

    net = DSGRN.Network(network)
    pg = DSGRN.ParameterGraph(net)

    one_change_hexes = {}
    for s,v in summary.items():
        try:
            w = int(v[0])
        except:
            continue
        if w > 0:
            hlist=[]
            for p in params[s]:
                param = pg.parameter(p)
                hexlist = [q.hex() for q in param.logic()]
                if sum([1 for h,d in zip(hexlist,design_spec) if h != d]) == 1:
                    hlist.append(hexlist)
            one_change_hexes.update({s : hlist})
    return summary,params,one_change_hexes


