import sys
sys.path.append('/home/jupyter/tacc-work/jupyter_packages/lib/python3.6/site-packages')
import DSGRN, graphviz
import progressbar, json, subprocess

def is_FP(annotation):
    return annotation.startswith("FP")

def is_FP_match(state, annotation):
    digits = [int(s) for s in annotation.replace(",", "").split() if s.isdigit()]
    return all(digits[k] >= state[k][0] and digits[k] <= state[k][1]
           for k in state)

def get_matching_truthtables(parametergraph,truthtables,N):
    params = [[] for _ in range(len(truthtables))]
    bar = progressbar.ProgressBar(max_value=N)
    for p in range(N):
        bar.update(p)
        parameter = parametergraph.parameter(p)
        dg = DSGRN.DomainGraph(parameter)
        md = DSGRN.MorseDecomposition(dg.digraph())
        mg = DSGRN.MorseGraph(dg, md)
        stable_FP_annotations = [mg.annotation(i)[0] for i in range(0, mg.poset().size())
                                 if is_FP(mg.annotation(i)[0]) and len(mg.poset().children(i)) == 0]
        for k,states in enumerate(truthtables):
            if all(any([is_FP_match(v,a) for a in stable_FP_annotations]) for v in states):
                params[k].append(p)
    bar.finish()
    return params

def results(net_str,truthtables,displaygraph=False):
    datetime = subprocess.check_output(['date +%Y_%m_%d_%H_%M_%S'], shell=True).decode(sys.stdout.encoding).strip()
    # truth tables is list of dictionaries
    network = DSGRN.Network()
    network.assign(net_str)
    if displaygraph:
        display(graphviz.Source(network.graphviz()))
    pg = DSGRN.ParameterGraph(network)
    np = pg.size()
    tt = []
    for states in truthtables:
        state_dicts = []
        for name, state in states.items():
            state_dicts.append({network.index(str(k)): state[k] for k in state})
        tt.append(tuple(state_dicts))
    params = get_matching_truthtables(pg,tt,np)
    for k,t in enumerate(tt):
        l = len(params[k])
        print("Truth table:")
        print(t)
        print("Parameters with specified truth table: {}/{} = {:.2f}%".format(l,np,l/np*100))
    D = {"network" : net_str}
    D.update( { k : t for k,t in enumerate(truthtables)} )
    json.dump(D, open('metadata{}.json'.format(datetime),'w'))
    json.dump(params, open('params{}.json'.format(datetime),'w'))
