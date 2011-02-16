import sqlite3
import re
import sys

# check arguments
if len(sys.argv) != 3:
	sys.exit("USAGE: python generateCharmmDB.py CharmTOP CharmPRM")

charmmtop=sys.argv[1]
charmmprm=sys.argv[2]

# OPEN NEW SQLITE3 DATABASE 'CHARMM.DB'
# Create table with the residue, atom name and atom type
# charge, radii and epsilon 
conn = sqlite3.connect('charmm.db')
c = conn.cursor()
c.execute("""
CREATE TABLE IF NOT EXISTS residues (
res_id INT AUTO_INCREMENT PRIMARY KEY,
name VARCHAR(10),
atname VARCHAR(4),
attype VARCHAR(4),
charge FLOAT,
radii FLOAT,
eps FLOAT);
""")
c.execute("CREATE INDEX name ON residues (name);")
c.execute("CREATE INDEX atname ON residues(atname)")
c.execute("CREATE INDEX attype ON residues (attype)")
conn.commit()

def insertResCharges(res, at):
	"Inserts into the database a new record"
        t = tuple([res]+at)
        c.execute("""INSERT INTO residues (name, atname, attype, charge) 
                   VALUES (?, ?, ?, ?)""",t)
        conn.commit()

# Open Topology Charmm file and parse information to the DB
f = open(charmmtop,'r')
l = True
atmatch=re.compile('^ATOM')
while(l):
        l=f.readline()
        if 'RESI' in l or 'PRES' in l:
                ats = []
                res = l.split()[1]
                print "inserting ",res
                while bool(l.strip()):
                        l = f.readline()
                        if atmatch.match(l):
                                s = l.split()
                                ats.append([s[1],s[2],s[3]])
                [insertResCharges(res, at) for at in ats]

f.close()

# Open Parameters Charmm file and add information to the previous records 
f=open(charmmprm,'r')

l=True
check=re.compile('^!')
t = []
while(l):
        l=f.readline()
	# Loop until finding !carbos (just a line marker for the section we want)
        if '!carbons' in l:
		# Loop and get information only from the lines not starting with !
                while (l):                        
			if not check.match(l):
                                s=l.split()
                                t.append([s[3],s[2],s[0]])                          
                        l = f.readline().strip()
f.close()

# Update information
c.executemany("UPDATE residues SET radii=?, eps=? WHERE attype=?",t)
conn.commit()

conn.close()
print "DONE"
