# The MIT License (MIT)
#
# Copyright (c) 2018 Bree Cummins
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import DSGRN,itertools,json,os,contextlib

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
    '''
    Makes the truth table look nice.
    '''
    l = []
    for t in tt:
        i1 = str(int(t[in1][0] > 0))
        i2 = str(int(t[in2][0] > 0))
        o = str(int(t[out][0] > 0))
        l.append(i1 + i2 + " " + o + "\n")
    return "".join(l)


def do_all_queries(dbname, in1, topval1, in2, topval2, out, topvalout, print_output=False):
    '''
    Take a database and record the parameters for each of the possible 16 truth tables for two inputs.

    :param dbname: String with the file name of a precomputed database via mpiexec Signatures.
    :param in1: String with the name of the first input node.
    :param topval1: The highest state (integer) that the input node can achieve.
    :param in2: Name of second input.
    :param topval2: Highest value of the second input.
    :param out: String with the name of the output variable.
    :param topvalout: Highest value of the output variable.
    :param print_output: True or False. Print the output for each truth table in pretty format

    :return: File names where the results are deposited. The keys to each json dictionary are
             integers that link truth tables in one file to parameters in the second file.
    '''
    database = DSGRN.Database(dbname)
    names = [database.names[k] for k in range(database.D)]
    for n in [in1,in2,out]:
        if n not in names:
            raise ValueError("{} is not a variable in the network for {}".format(n,dbname))
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
    sum_file = 'summary_{}.json'.format(dbname[:-3])
    param_file = 'params_{}.json'.format(dbname[:-3])
    D = {"database" : dbname, "num_params" : np }
    D.update( { k : (l,t) for (k,t),l in zip(enumerate(truthtables),lengths) } )
    json.dump(D, open(sum_file,'w'))
    json.dump({ k : a for k,a in enumerate(all_matches)}, open(param_file,'w'))
    return sum_file,param_file
