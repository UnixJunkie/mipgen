#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="dalvarez"
__date__ ="$17-feb-2011 22:52:41$"

def insertAt(type, r, e):
    t = (type, r, e)
    c.execute("INSERT INTO atypes VALUES (null, ?, ?, ?)",t)
    conn.commit()

def parseMOD4(lines):
    mod4 = re.compile('^MOD4')
    for i,l in enumerate(lines):
        if mod4.match(l):
            N_dict = lines[i-3]
            C_dict = lines[i-2]
            print N_dict
            print C_dict
            c = 1
            l = lines[i+c].split()
            while(bool(l)):
                insertAt(l[0],l[1],l[2])
                c += 1
                l = lines[i+c].split()
    return N_dict, C_dict

if __name__ == "__main__":
    import sqlite3
    import sys
    import re

    if len(sys.argv) != 3:
        sys.exit("\n\tUSAGE: python generateAMBERDB.py parmXX.dat gaff.dat")

    dat = sys.argv[1]
    gaff = sys.argv[2]

    # OPEN NEW SQLITE3 DATABASE 'AMBER.DB'
    # Create table with atom type, and VdW radii and epsilon
    conn = sqlite3.connect('amber.db')
    c = conn.cursor()
    c.execute("""
    CREATE TABLE IF NOT EXISTS atypes (
    at_id INT AUTO_INCREMENT PRIMARY KEY,
    attype VARCHAR(4),
    radii FLOAT,
    eps FLOAT);
    """)
    c.execute("CREATE INDEX attype ON atypes (attype)")
    conn.commit()

    # PARSE PARM.dat FILE
    # Find MOD4 section and read 2 previous lines also
    lines = [line.strip() for line in open(dat, 'r').readlines()]
    N_dict, C_dict = parseMOD4(lines)
    
    # Add atoms in the N_dict and C_dict list also
    # with the parameters of the first element which where
    # added in the previous step
    N_dict = N_dict.split()
    C_dict = C_dict.split()
    for d in (N_dict, C_dict):
        at = d.pop(0)
        for a in d:
            c.execute("INSERT INTO atypes(attype) VALUES('%s')"%a)
            c.execute("UPDATE atypes SET radii = (SELECT radii FROM atypes WHERE attype=?) WHERE attype=?",(at, a))
            c.execute("UPDATE atypes SET eps = (SELECT eps FROM atypes WHERE attype=?) WHERE attype=?",(at, a))
            conn.commit()

    # PARSE GAFF.dat FILE
    # Same way as before but without C and N dict list
    # so the output will be ignored
    lines = [line.strip() for line in open(gaff, 'r').readlines()]
    N_dict, C_dict = parseMOD4(lines)

    # Close connections to the DB
    # and we are done!
    conn.close()
    print "DONE"