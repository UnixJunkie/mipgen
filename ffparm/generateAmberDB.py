#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.
__author__="dalvarez"
__date__ ="$17-feb-2011 22:52:41$"

def insertAt(*args):
    "args is: (resname, atname, attype, charge) in this strict order"
    c.execute("INSERT INTO ambertypes(resname, atname, attype, charge) VALUES (?, ?, ?, ?)",args)
    conn.commit()

def insertVdW(*args):
    "insert into the table the VdW parameters for a given atomtype.\
    *args: (attype, radii, eps)"
    c.execute("UPDATE ambertypes SET radii=%.4f, eps=%.4f WHERE attype='%s'"%(args[1],args[2],args[0]))
    conn.commit()

def parseMOD4(dat):
    lines = [line.strip() for line in open(dat, 'r').readlines()]
    mod4 = re.compile('^MOD4')
    for i,l in enumerate(lines):
        if mod4.match(l):
#            N_dict = lines[i-3]
#            C_dict = lines[i-2]
            c = 1
            l = lines[i+c].split()
            while(bool(l)):
                attype = l[0]
#                if attype == el for el in ['N','NA','N2','N*','NC','NB','NT','NY']
                insertVdW(l[0],l[1],l[2])
                c += 1
                l = lines[i+c].split()
                
#    return N_dict, C_dict

def parseLibFile(file):
    """Same format as a Prepin file"""
    from amberprep import readAmberPrep
    res_dict = readAmberPrep(file)
    # Add residues, atom names, atomtypes and charges to the DB
    [insertAt(res,at[0],at[1],at[2]) for at in res_dict[res] for res in res_dict.keys()]
    return True

if __name__ == "__main__":
    import sqlite3
    import sys
    import re

    if len(sys.argv) != 3:
        sys.exit("\n\tUSAGE: python generateAMBERDB.py all_aminoXX.in parmXX.dat ")

    dat = sys.argv[2]
    lib = sys.argv[1]

    # OPEN NEW SQLITE3 DATABASE 'AMBER.DB'
    # Create table with atom type, and VdW radii and epsilon
    conn = sqlite3.connect('amber.db')
    c = conn.cursor()
    c.execute("""
    CREATE TABLE IF NOT EXISTS ambertypes (
    at_id INT AUTO_INCREMENT PRIMARY KEY,
    resname VARCHAR(5),
    atname VARCHAR(5),
    attype VARCHAR(4),
    charge FLOAT,
    radii FLOAT,
    eps FLOAT);
    """)
    c.execute("CREATE INDEX attype ON ambertypes (attype)")
    c.execute("CREATE INDEX atname ON ambertypes (atname)")
    c.execute("CREATE INDEX resname ON ambertypes (resname)")
    conn.commit()

    # PARSE ALL_AMINOXX.IN (same format as a Prepin File)
    # This file contains all the residue names supported by amberFF
    # with a atom name to atom type mapping and the charges for each atomtype
    parseLibFile(lib)
    
    # PARSE PARM.dat FILE
    # Find MOD4 section and read 2 previous lines also
#    N_dict, C_dict = parseMOD4(dat)
    parseMOD4(dat)
    
    # Add atoms in the N_dict and C_dict list also
    # with the parameters of the first element which where
    # added in the previous step
#    N_dict = N_dict.split()
#    C_dict = C_dict.split()
#    for d in (N_dict, C_dict):
#        at = d.pop(0)
#        for a in d:
#            c.execute("INSERT INTO ambertypes(attype) VALUES('%s')"%a)
#            c.execute("UPDATE atypes SET radii = (SELECT radii FROM atypes WHERE attype=?) WHERE attype=?",(at, a))
#            c.execute("UPDATE atypes SET eps = (SELECT eps FROM atypes WHERE attype=?) WHERE attype=?",(at, a))
#            conn.commit()

    # MANUALLY ADD H0 Hydrogen for GLY
    # incorporated in ffrcmod.03 and used nowadays
    #  H0       1.3870   0.0157             Veenstra et al JCC,8,(1992),963
    insertAt('H0',1.387,0.0157)
    
    # Close connections to the DB
    # and we are done!
    conn.close()
    print "DONE"