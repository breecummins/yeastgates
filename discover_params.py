import DSGRN,json,sqlite3
from truth_table_fixed_pts import do_all_queries


def NAND():
    dbname = "NAND.db"
    topval1 = 2
    topval2 = 2
    design_spec = ["C","C","2","2","8","2"]
    in1 = "s1"
    in2 = "s2"
    return dbname, in1, topval1, in2, topval2, design_spec


def OR():
    dbname = "OR.db"
    topval1 = 2
    topval2 = 2
    design_spec = ["C","C","8","2"]
    in1 = "s1"
    in2 = "s2"
    return dbname, in1, topval1, in2, topval2, design_spec


def AND():
    dbname = "AND.db"
    topval1 = 2
    topval2 = 2
    design_spec = ["C","C","2","2","8"]
    in1 = "s1"
    in2 = "s2"
    return dbname, in1, topval1, in2, topval2, design_spec


def XOR():
    dbname = "XOR.db"
    topval1 = 3
    topval2 = 3
    design_spec = ["38","38","8","2","2","8","8"]
    in1 = "s1"
    in2 = "s2"
    return dbname, in1, topval1, in2, topval2, design_spec


def get_network(dbname):
    conn = sqlite3.connect(dbname)
    cursor = conn.cursor()
    cursor.execute("select Specification from Network")
    network = cursor.fetchone()[0]
    conn.close()
    return network


def get_hex(pind,pg):
    return [q.hex() for q in pg.parameter(pind).logic()]


def get_one_change_params(logic,out="GFP",topvalout=1,print_output=True):
    '''

    :param logic: A function handle specifying one of the logics in this module.
    :return: A dictionary of truth tables and parameters counts, a corresponding dictionary of the parameter indices, and a dictionary of the most parsimonious hex codes for each truth table (the hex codes where there is exactly one change from design spec.
    '''
    dbname, in1, topval1, in2, topval2, design_spec = logic()
    sfile, pfile = do_all_queries(dbname, in1, topval1, in2, topval2, out, topvalout, print_output)

    summary = json.load(open(sfile))
    params = json.load(open(pfile))

    pg = DSGRN.ParameterGraph(DSGRN.Network(get_network(dbname)))

    one_change_hexes = {}
    for s,v in summary.items():
        try:
            w = int(v[0])
        except:
            continue
        if w > 0:
            hlist=[]
            for p in params[s]:
                hexlist = get_hex(p,pg)
                if sum([1 for h,d in zip(hexlist,design_spec) if h != d]) == 1:
                    hlist.append(hexlist)
            one_change_hexes.update({s : hlist})
    return summary,params,one_change_hexes


