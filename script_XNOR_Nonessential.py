from yeastgatemodule import results

# XNOR gate
with open("XNOR.txt",'r') as f:
    net_str = f.read()

desiredtt = {
    "state00" : {"s1":[0,0],"s2":[0,0],"GFP":[1,1]},
    "state01" : {"s1":[0,0],"s2":[3,3],"GFP":[0,0]},
    "state10" : {"s1":[3,3],"s2":[0,0],"GFP":[0,0]},
    "state11" : {"s1":[3,3],"s2":[3,3],"GFP":[1,1]}
}

fail1 = {
    "state00" : {"s1":[0,0],"s2":[0,0],"GFP":[0,0]},
    "state01" : {"s1":[0,0],"s2":[3,3],"GFP":[1,1]},
    "state10" : {"s1":[3,3],"s2":[0,0],"GFP":[1,1]},
    "state11" : {"s1":[3,3],"s2":[3,3],"GFP":[1,1]}
}

fail2 = {
    "state00" : {"s1":[0,0],"s2":[0,0],"GFP":[0,0]},
    "state01" : {"s1":[0,0],"s2":[3,3],"GFP":[0,0]},
    "state10" : {"s1":[3,3],"s2":[0,0],"GFP":[1,1]},
    "state11" : {"s1":[3,3],"s2":[3,3],"GFP":[0,0]}
}

# should be the same as fail2
fail3 = {
    "state00" : {"s1":[0,0],"s2":[0,0],"GFP":[0,0]},
    "state01" : {"s1":[0,0],"s2":[3,3],"GFP":[1,1]},
    "state10" : {"s1":[3,3],"s2":[0,0],"GFP":[0,0]},
    "state11" : {"s1":[3,3],"s2":[3,3],"GFP":[0,0]}
}

fail4 = {
    "state00" : {"s1":[0,0],"s2":[0,0],"GFP":[1,1]},
    "state01" : {"s1":[0,0],"s2":[3,3],"GFP":[0,0]},
    "state10" : {"s1":[3,3],"s2":[0,0],"GFP":[0,0]},
    "state11" : {"s1":[3,3],"s2":[3,3],"GFP":[0,0]}
}

fail5 = {
    "state00" : {"s1":[0,0],"s2":[0,0],"GFP":[1,1]},
    "state01" : {"s1":[0,0],"s2":[3,3],"GFP":[1,1]},
    "state10" : {"s1":[3,3],"s2":[0,0],"GFP":[0,0]},
    "state11" : {"s1":[3,3],"s2":[3,3],"GFP":[0,0]}
}

fail6 = {
    "state00" : {"s1":[0,0],"s2":[0,0],"GFP":[1,1]},
    "state01" : {"s1":[0,0],"s2":[3,3],"GFP":[0,0]},
    "state10" : {"s1":[3,3],"s2":[0,0],"GFP":[1,1]},
    "state11" : {"s1":[3,3],"s2":[3,3],"GFP":[0,0]}
}

fail7 = {
    "state00" : {"s1":[0,0],"s2":[0,0],"GFP":[0,0]},
    "state01" : {"s1":[0,0],"s2":[3,3],"GFP":[1,1]},
    "state10" : {"s1":[3,3],"s2":[0,0],"GFP":[1,1]},
    "state11" : {"s1":[3,3],"s2":[3,3],"GFP":[0,0]}
}

fail8 = {
    "state00" : {"s1":[0,0],"s2":[0,0],"GFP":[0,0]},
    "state01" : {"s1":[0,0],"s2":[3,3],"GFP":[1,1]},
    "state10" : {"s1":[3,3],"s2":[0,0],"GFP":[0,0]},
    "state11" : {"s1":[3,3],"s2":[3,3],"GFP":[1,1]}
}

fail9 = {
    "state00" : {"s1":[0,0],"s2":[0,0],"GFP":[0,0]},
    "state01" : {"s1":[0,0],"s2":[3,3],"GFP":[0,0]},
    "state10" : {"s1":[3,3],"s2":[0,0],"GFP":[1,1]},
    "state11" : {"s1":[3,3],"s2":[3,3],"GFP":[1,1]}
}

fail10 = {
    "state00" : {"s1":[0,0],"s2":[0,0],"GFP":[1,1]},
    "state01" : {"s1":[0,0],"s2":[3,3],"GFP":[1,1]},
    "state10" : {"s1":[3,3],"s2":[0,0],"GFP":[1,1]},
    "state11" : {"s1":[3,3],"s2":[3,3],"GFP":[0,0]}
}

fail11 = {
    "state00" : {"s1":[0,0],"s2":[0,0],"GFP":[1,1]},
    "state01" : {"s1":[0,0],"s2":[3,3],"GFP":[1,1]},
    "state10" : {"s1":[3,3],"s2":[0,0],"GFP":[0,0]},
    "state11" : {"s1":[3,3],"s2":[3,3],"GFP":[1,1]}
}

fail12 = {
    "state00" : {"s1":[0,0],"s2":[0,0],"GFP":[0,0]},
    "state01" : {"s1":[0,0],"s2":[3,3],"GFP":[0,0]},
    "state10" : {"s1":[3,3],"s2":[0,0],"GFP":[0,0]},
    "state11" : {"s1":[3,3],"s2":[3,3],"GFP":[0,0]}
}

fail13 = {
    "state00" : {"s1":[0,0],"s2":[0,0],"GFP":[1,1]},
    "state01" : {"s1":[0,0],"s2":[3,3],"GFP":[1,1]},
    "state10" : {"s1":[3,3],"s2":[0,0],"GFP":[1,1]},
    "state11" : {"s1":[3,3],"s2":[3,3],"GFP":[1,1]}
}

fail14 = {
    "state00" : {"s1":[0,0],"s2":[0,0],"GFP":[0,0]},
    "state01" : {"s1":[0,0],"s2":[3,3],"GFP":[0,0]},
    "state10" : {"s1":[3,3],"s2":[0,0],"GFP":[0,0]},
    "state11" : {"s1":[3,3],"s2":[3,3],"GFP":[1,1]}
}

fail15 = {
    "state00" : {"s1":[0,0],"s2":[0,0],"GFP":[1,1]},
    "state01" : {"s1":[0,0],"s2":[3,3],"GFP":[0,0]},
    "state10" : {"s1":[3,3],"s2":[0,0],"GFP":[1,1]},
    "state11" : {"s1":[3,3],"s2":[3,3],"GFP":[1,1]}
}

results(net_str,"XNOR",[desiredtt, fail1, fail2, fail3, fail4, fail5, fail6, fail7, fail8, fail9, fail10, fail11, fail12, fail13, fail14, fail15])