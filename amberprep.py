# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="dalvarez"
__date__ ="$18-feb-2011 0:14:26$"

import re

def readAmberPrep(prepfile):
    "Read an Amber prepin file and return the residues, atom names, atom types,\
    and charges for each atom in a dictionary {RESNAME:[[atname1,attype1,charge1], [..],..], ...}"
    res = re.compile(r'^\w+.res')       # Begin of a residue in the prepin file
    molec = {}

    f = open(prepfile, 'r')
    l = f.readline()
    while(l):
        l = f.readline()
        if res.match(l):
            resname = f.readline().split()[0]
            molec[resname] = []
            #skip 2 lines
            f.readline()
            f.readline()

            # begins atom section
            s = f.readline().strip()
            while(bool(s)):
                atline = s.split()
                atname, attype, charge = atline[1], atline[2], atline[-1]
                if atname != 'DUMM': molec[resname].append((atname, attype, charge))
                s = f.readline().strip()
                
    f.close()
    return molec

if __name__ == "__main__":
    import sys
    molec = readAmberPrep(sys.argv[1])
    for res in molec.keys():
        print res, molec[res]
