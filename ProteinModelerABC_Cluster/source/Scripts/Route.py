import argparse

parser = argparse.ArgumentParser()

parser.add_argument("route")

arg = parser.parse_args()

# Installation route
def route_bin(route):

    route = "'" + route + "'"
    
    with open('ProteinModelerABC_Cluster.py', 'r') as f:
        lines = f.readlines()

    with open('ProteinModelerABC_Cluster.py', 'w') as nf:
        for line in lines:
            if line.startswith('sys.path.append('):
                nf.write('sys.path.append(' + str(route) + ')\n')
            elif line.startswith('ROOT = '):
                nf.write('ROOT = "' + str(route) + '"\n')
            else:
                nf.write(line)
                
                
    with open('./bin/Scripts/Functions.py', 'r') as g:
        ls = g.readlines()
    
    with open('./bin/Scripts/Functions.py', 'w') as ng:
        for l in ls:
            if 'INSTALLATIONROOT' in l:
                l = l.replace('INSTALLATIONROOT', route)
            ng.write(l)            

                
route_bin(arg.route)