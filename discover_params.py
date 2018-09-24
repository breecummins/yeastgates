import DSGRN,json,sqlite3
from truth_table_fixed_pts import do_all_queries


def dict_of_behavior(node_type,hex):
    # node type = (num in-edges, num out-edges)
    # behaviors = ["Design","Cellular repression","Cellular activation","Leaky repression at 1 promoter","Leaky repression at 2 promoters","Mixed effects"]

    if hex == 0:
        return "Cellular repression"

    if node_type == (1,1):
        if hex == 2:
            return "Design"
        elif hex == 3:
            return "Cellular activation or leaky repression"
    elif node_type == (2,1):
        if hex == 8:
            return "Design"
        elif hex == 0x0A or 0x0C:
            return "Cellular activation or leaky repression at one promoter"
        elif hex == 0x0E or 0x0F:
            return "Cellular activation or leaky repression at two promoters"
    elif node_type == (1,2):
        #This works only for an ACTIVATING in-edge, as for either of the input nodes to the circuit
        if hex == 4:
            return "Cellular repression"
        elif hex == 5:
            return "Mixed effects"
        elif hex == 0x0C:
            return "Design"
        elif hex == 0x0D or 0x0F:
            return "Cellular activation"
    elif node_type not in [(1,1),(2,1),(1,2)]:
        return "Node type {} is not implemented.".format(node_type)

    return "{} is not a hex code for a {} node type.".format(hex,node_type)



def NOR():
    dbname = "NOR.db"
    topval1 = 2
    topval2 = 2
    node_types = [(1,2),(1,2),(2,1)]
    design_spec = ["C","C","8"]
    in1 = "s1"
    in2 = "s2"
    return dbname, in1, topval1, in2, topval2, design_spec, node_types


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


def XNOR():
    #FIXME: Design spec is incomplete
    dbname = "XNOR.db"
    topval1 = 3
    topval2 = 3
    design_spec = ["38","38","","8","8","8"]
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


