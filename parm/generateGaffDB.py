#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="dalvarez"
__date__ ="$21-feb-2011 19:17:23$"

def insertAt(type, r, e):
    t = (type, r, e)
    c.execute("INSERT INTO gafftypes VALUES (null, ?, ?, ?)",t)
    conn.commit()

def parseMOD4(dat):
    lines = [line.strip() for line in open(dat, 'r').readlines()]
    mod4 = re.compile('^MOD4')
    for i,l in enumerate(lines):
        if mod4.match(l):
            N_dict = lines[i-3]
            C_dict = lines[i-2]
            c = 1
            l = lines[i+c].split()
            while(bool(l)):
                attype = l[0]
                insertAt(attype,l[1],l[2])
                c += 1
                l = lines[i+c].split()

#    return N_dict, C_dict

if __name__ == "__main__":
    import sqlite3
    import sys
    import re

    if len(sys.argv) != 2:
        sys.exit("\n\tUSAGE: python generateGaffDB.py gaff.dat")

    gaff = sys.argv[1]

    # OPEN NEW SQLITE3 DATABASE 'AMBER.DB'
    # Create table with atom type, and VdW radii and epsilon
    conn = sqlite3.connect('gaff.db')
    c = conn.cursor()
    c.execute("""
    CREATE TABLE IF NOT EXISTS gafftypes (
    at_id INT AUTO_INCREMENT PRIMARY KEY,
    attype VARCHAR(4),
    radii FLOAT,
    eps FLOAT);
    """)
    c.execute("CREATE INDEX attype ON gafftypes (attype)")
    conn.commit()

    # PARSE GAFF.dat FILE
    # Same way as before but without C and N dict list
    # so the output will be ignored
    parseMOD4(gaff)

    # Close connections to the DB
    # and we are done!
    conn.close()
    print "DONE"