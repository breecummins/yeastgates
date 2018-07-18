import sys,os

fname = sys.argv[1]

with open(fname,'r') as f:
    f.readline()
    newname = os.path.splitext(fname)[0] + "_links.txt"
    links = []
    for l in f:
        words = l.split(',')
        if words[3] == "Audit Passed":
            links.append(words[4])
    with open(newname,'w') as n:
        line = '["' + '","'.join(links) + '"]'
        n.write(line)

