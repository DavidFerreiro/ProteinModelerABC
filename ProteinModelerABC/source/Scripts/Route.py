import argparse

parser = argparse.ArgumentParser()

parser.add_argument("route")

arg = parser.parse_args()

# Installation route
def route_bin(route):
    with open('ProteinModelerABC.py', 'r') as f:
        lines = f.readlines()

    with open('ProteinModelerABC.py', 'w') as nf:
        for line in lines:
            if line.startswith('sys.path.append('):
                nf.write('sys.path.append("' + str(route) + '")\n')
            else:
                nf.write(line)
                
    route = "'" + route + "'"
    with open('./bin/Scripts/Variables.py', 'r') as g:
        ls = g.readlines()
    
    with open('./bin/Scripts/Variables.py', 'w') as ng:
        for l in ls:
            if l.startswith('INSTALLATIONROOT ='):
                ng.write('INSTALLATIONROOT =' + str(route) + '\n')   
            else:
                ng.write(l)       

                
route_bin(arg.route)