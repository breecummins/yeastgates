from yeastgatemodule import results

# XOR gate
# This topology makes up for the fact that we don't have AB + CD logic in DSGRN.
net_str = "s1 : s1\ns2 : s2\ns3 : (~s1)(~s2)\ns4 : ~s2\ns5 : ~s1\ns6 : (~s4)(~s5)\nGFP : ~s3 + ~s6\nD : (D)(GFP) : E"

desiredtt = {
    "state00" : {"s1":[0,0],"s2":[0,0],"GFP":[0,0]},
    "state01" : {"s1":[0,0],"s2":[1,3],"GFP":[1,1]},
    "state10" : {"s1":[1,3],"s2":[0,0],"GFP":[1,1]},
    "state11" : {"s1":[1,3],"s2":[1,3],"GFP":[0,0]}
}

fail1 = {
    "state00" : {"s1":[0,0],"s2":[0,0],"GFP":[0,0]},
    "state01" : {"s1":[0,0],"s2":[1,3],"GFP":[1,1]},
    "state10" : {"s1":[1,3],"s2":[0,0],"GFP":[1,1]},
    "state11" : {"s1":[1,3],"s2":[1,3],"GFP":[1,1]}
}

fail2 = {
    "state00" : {"s1":[0,0],"s2":[0,0],"GFP":[0,0]},
    "state01" : {"s1":[0,0],"s2":[1,3],"GFP":[0,0]},
    "state10" : {"s1":[1,3],"s2":[0,0],"GFP":[1,1]},
    "state11" : {"s1":[1,3],"s2":[1,3],"GFP":[0,0]}
}

# should be the same as fail2
fail3 = {
    "state00" : {"s1":[0,0],"s2":[0,0],"GFP":[0,0]},
    "state01" : {"s1":[0,0],"s2":[1,3],"GFP":[1,1]},
    "state10" : {"s1":[1,3],"s2":[0,0],"GFP":[0,0]},
    "state11" : {"s1":[1,3],"s2":[1,3],"GFP":[0,0]}
}
results(net_str,[desiredtt, fail1, fail2, fail3])