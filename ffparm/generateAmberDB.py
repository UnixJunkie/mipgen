#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.
__author__="dalvarez"
__date__ ="$17-feb-2011 22:52:41$"

import re

def insertAt(*args):
    "args is: (resname, atname, attype, charge) in this strict order"
    args = tuple(args)
    c.execute("INSERT INTO ambertypes(resname, atname, attype, charge) VALUES (?, ?, ?, ?)",args)
    conn.commit()

def insertVdW(*args):
    "insert into the table the VdW parameters for a given atomtype.\
    *args: (attype, radii, eps)"
    c.execute("UPDATE ambertypes SET radii=%s, eps=%s WHERE attype='%s'"%(args[1],args[2],args[0]))
    conn.commit()

def parseAmberParm(dat):
    "Reads an Amber Parm file and extracts the MOD4 section. The one containing\
    VdW parameters for each atomtype. 'dat' is the name of the file (string)."
    lines = [line.strip() for line in open(dat, 'r').readlines()]

    # This indicates the beggining of the section
    mod4 = re.compile('^MOD4')

    # Loop over the file lines and stop at MOD4
    for i,l in enumerate(lines):
        if mod4.match(l):

            # Take also 2 previous lines containing mapping names
            # These N and C atoms map all to the same parameters
            N_dict = lines[i-3].split()
            C_dict = lines[i-2].split()
            c = 1
            l = lines[i+c].split()
            while(bool(l)):
                attype = l[0]
                if attype == 'N':
                    # Assign same parameters for all of the names in N_dict
                    for el in N_dict:
                        insertVdW(el,l[1],l[2])
                elif attype == 'C*':
                    # Same for C_dict
                    for el in C_dict:
                        insertVdW(el,l[1],l[2])
                else:
                    # Other atoms
                    insertVdW(attype,l[1],l[2])
                c += 1
                l = lines[i+c].split()
                
    return True

def parseLibFile(file):
    """Same format as a Prepin file"""
    from amberprep import readAmberPrep
    res_dict = readAmberPrep(file)
    # Add residues, atom names, atomtypes and charges to the DB
    for res in res_dict.keys():
        print "Parsing: ",res
        [insertAt(res,at[0],at[1],at[2]) for at in res_dict[res]]
        
    return True

if __name__ == "__main__":
    import sqlite3
    import sys
    import re
    
    if len(sys.argv) != 5:
        sys.exit("\n\tUSAGE: python generateAMBERDB.py all_aminoXX.in all_aminoctXX.in all_aminontXX.in parmXX.dat\n\n")
        
    all = sys.argv[1]
    all_ct = sys.argv[2]
    all_nt = sys.argv[3]
    dat = sys.argv[4]
    
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
    
    # PARSE ALL_AMINOXXX.IN (same format as a Prepin File)
    # These file contains all the residue names supported by amberFF
    # with a atom name to atom type mapping and the charges for each atomtype

    # 3 files are needed for a FF version: all_amino03.in (standard aa params)
    # all_aminoct03.in (C-terminus aa) and all_aminont03.in (N-terminus aa)
    parseLibFile(all)
    parseLibFile(all_ct)
    parseLibFile(all_nt)
    
    # PARSE PARM.dat FILE
    # Find MOD4 section and read VdW parameters
    parseAmberParm(dat)

    # MANUALLY ADD H0 Hydrogen for GLY
    # incorporated in ffrcmod.03 and used nowadays
    #  H0       1.3870   0.0157             Veenstra et al JCC,8,(1992),963
    insertVdW('H0',1.387,0.0157)
    
    # Close connections to the DB
    # and we are done!
    conn.close()
    print "DONE"