import DSGRN,itertools,subprocess,sys,json,os,contextlib

def truth_table_constructor_2in(in1, topval1, in2, topval2, out, topvalout):
    '''

    :param in1: string with network var representing input 1
    :param topval1: top state in network for input 1
    :param in2: analogous for input 2
    :param topval2: analogous for input 2
    :param out: analogous for output
    :param topvalout: analogous for output
    :return: list of lists of truth table dictionaries
    '''
    truthtables = []
    combs = [c  for k in range(5) for c in itertools.combinations(range(4), k)]
    for c in combs:
        temp = [
        {in1: [0, 0], in2: [0, 0], out: [0, 0]},
        {in1: [0, 0], in2: [topval2,topval2], out: [0, 0]},
        {in1: [topval1, topval1], in2: [0, 0], out: [0, 0]},
        {in1: [topval1, topval1], in2: [topval2, topval2], out: [0, 0]}
        ]
        for ind in c:
            temp[ind][out] = [topvalout, topvalout]
        truthtables.append(tuple(temp))
    return truthtables

def truth_table_query(database,bounds):
    '''
    This function only works for fixed points with non-overlapping bounds.
    See DSGRN.DoubleFixedPointQuery for a technique when there are overlapping bounds.
    :param database: A DSGRN database.
    :param bounds: A list of bounds in DSGRN dictionary format.
    :return: a list of parameter indices with the matching truth table
    '''
    #FIXME: Add check that bounds do not overlap (they shouldn't for truth tables)
    for k,bound in enumerate(bounds):
        with open(os.devnull, 'w') as devnull:
            with contextlib.redirect_stdout(devnull):
                DSGRN.Query.FixedPointTables.MatchQuery(bound, "tempFP", database)
        c = database.conn.cursor()
        sqlexpression = 'select ParameterIndex from tempFP natural join Signatures;'
        paramlist = set([row[0] for row in c.execute(sqlexpression)])
        if not k:
            matches = paramlist
        else:
            matches = matches.intersection(paramlist)
        c.execute("drop table tempFP")
    return matches


def format_truth_table(tt,in1,in2,out):
    l = []
    for t in tt:
        i1 = str(int(t[in1][0] > 0))
        i2 = str(int(t[in2][0] > 0))
        o = str(int(t[out][0] > 0))
        l.append(i1 + i2 + " " + o + "\n")
    return "".join(l)


def do_all_queries(dbname, in1, topval1, in2, topval2, out, topvalout, print_output=False):
    database = DSGRN.Database(dbname)
    truthtables = truth_table_constructor_2in(in1, topval1, in2, topval2, out, topvalout)
    all_matches = []
    lengths = []
    for bounds in truthtables:
        matches = truth_table_query(database, bounds)
        all_matches.append(sorted(list(matches)))
        lengths.append(len(matches))
    np = database.parametergraph.size()
    if print_output:
        for l, tt in sorted(zip(lengths,truthtables),key=lambda t : t[0],reverse=True):
            if l > 0:
                print("Parameters with specified truth table: {}/{} = {:.2f}%\n".format(l, np,l/np*100))
                print(format_truth_table(tt,in1,in2,out))
    datetime = subprocess.check_output(['date +%Y_%m_%d_%H_%M_%S'], shell=True).decode(sys.stdout.encoding).strip()
    D = {"database" : dbname, "num_params" : np }
    D.update( { k : (l,t) for (k,t),l in zip(enumerate(truthtables),lengths) } )
    json.dump(D, open('summary{}_{}.json'.format(datetime,dbname[:-3]),'w'))
    json.dump({ k : a for k,a in enumerate(all_matches)}, open('params{}_{}.json'.format(datetime,dbname[:-3]),'w'))



# def truth_table_query_old(database, bounds):
#     '''
#     This class only works for fixed points with non-overlapping bounds.
#     See DSGRN.DoubleFixedPointQuery for a technique when there are overlapping bounds.
#     :param database: A DSGRN database.
#     :param bounds: A list of bounds in DSGRN dictionary format.
#     :return: a list of parameter indices with the matching truth table
#     '''
#     #FIXME: Add check that bounds do not overlap (they shouldn't for truth tables)
#     tables = ["Matches{}".format(k) for k in range(len(bounds))]
#     for name, b in zip(tables, bounds):
#         DSGRN.Query.FixedPointTables.MatchQuery(b, name, database);
#     make_query = "".join([" intersect select * from {}".format(n) for n in tables[1:]])
#     c = database.conn.cursor()
#     c.execute("create temp table Matches as select MorseGraphIndex from (select * from {}".format(tables[0]) + make_query + ") group by MorseGraphIndex;")
#     matches = [row[0] for row in c.execute('select distinct ParameterIndex from Matches natural join Signatures;')]
#     c.execute('drop table Matches')
#     for name in tables:
#         c.execute('drop table {}'.format(name))
#     return matches


