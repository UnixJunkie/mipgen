#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="dalvarez"
__date__ ="$23-feb-2011 18:22:26$"
import re

def readResAtCoordFromPDB(pdbFile):
    atom = re.compile(r'^(ATOM|HETATM)\s')
    f = open(pdbFile, 'r')
    l = f.readline()
    pdb = []
    while l:
        if atom.match(l):
            # ATOMNAME //13-16 columns
            at = l[12:16].strip()
            # If atom name begins with a number, put number at the end
            if re.match(r'\d',at[0]): at = at[1:]+at[0] 

            # RESIDUE NAME // 18-20 cols
            res = l[17:20].strip()
            
            # COORDINATES // 31-54 cols
            xyz = map(float, l[30:54].split())
            pdb.append([res,at]+xyz)
            
        l = f.readline()
    f.close()
    return pdb

if __name__ == "__main__":
    import sys
    print readResAtCoordFromPDB(sys.argv[1])
