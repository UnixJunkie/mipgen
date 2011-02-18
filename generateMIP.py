#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="dalvarez"
__date__ ="$04-feb-2011 18:39:53$"

import Biskit as bi
import sqlite3
import numpy as npy

# Setting up some constants
vdw_cutoff=10.
elec_cutoff=20.
dielectric_constant = 80.
coul_k = 332./(dielectric_constant)

# Fake PROBE parameters for CHARGE,VDWradii,Epsilon
probeparms = {'C':[0.0,1.95,-0.07], # CD from leucine
             'O':[-0.01,1.65,-0.12],
             'H+':[1,1.1,-0.0157],
             'N':[-0.47,1.85,-0.2],
             'S':[-0.09,2.1,-0.47],
             'neg':[-1,2.5,-0.7]}

def distanceMat(xyz1, xyz2=None):
    """Calculate distance matrix.
    If only xyz1 set is given. Calculate DM with itself.
    If xyz2 is give, it should be only a point. Then calculate DM between the point and xyz1.
    """
    
    if xyz2 is not None:
        n1 = len(xyz1)
        dm = npy.zeros((1, n1))
        dm = ((xyz1 - xyz2)**2).sum(axis=1)
        
    else:
        n = len(xyz1)
        dm = npy.zeros((n,n))
        for i in range(3):
            data = xyz1[:,i]
            dm += (data - data[:,npy.newaxis])**2
    return npy.sqrt(dm)


class Molecule:

    def __init__(self, pdbModel):
        """Build an instance for the molecule containing
        for each atom, the coordinates and non-bonded parameters
        (mass, charge, vdw radii and epsilon)

        Input:  atom_names - string with names as in PDB
                xyz - numpy.ndarray with the coordinates
        """
        
        self.atoms, self.missing = self.__getParamsXYZ(pdbModel)
        self.natoms = len(self.atoms)
        self.xyz = npy.array([at[0][3:] for at in self.atoms])
            

    def __getParamsXYZ(self, pdbInstance):
        """Returns for each atom in pdbInstance a tuple (charge, radii, epsilon, x, y, z)"""

        # Some name conversions
        resDict={'WAT':'TIP3','HOH':'TIP3','HIS':'HSE','HIE':'HSE','HID':'HSD','HIP':'HSP'}
        atlist={'ILE':{'CD1':'CD'}}
        
        # Open sqlite3 database connection to charmm DB
        conn = sqlite3.connect('ffparm/charmm.db')
        c = conn.cursor()

        # Initiate empty lists
        # missing list will contain the atom indices not recognized in the DB
        atoms = []
        missing = []
        
        # sotre coordinates, residue names and atom names 
        # for faster access during the loop
        xyz = pdbInstance.xyz
        resnames = pdbInstance['residue_name']
        atnames = pdbInstance['name']
        
        for i in range(len(pdbInstance)):
            res,at = resnames[i], atnames[i]
            if res in resDict.keys(): res=resDict[res]
            if res in atlist.keys():
                if at in atlist[res].keys(): at=atlist[res][at]

            # Fetch parameters
            c.execute("SELECT charge,radii,eps FROM residues WHERE name=? AND atname=?",(res,at))
            row=c.fetchall()
            if len(row) != 1:
                # No results. Append to missing list
                print "Not found: Atom %i res %s atom %s. Or multiple results"%(i,res,at)
                missing.append(i)
            else:
                # Correct fetching
                row=map(float, row[0])
                atoms.append([[row[0],row[1],row[2]]+xyz[i].tolist()])
        conn.close()

        return atoms, missing

    def __getitem__(self, idx):
        return self.atoms[idx][0]


def calcVdW(probe, mol, dmat):
    V = 0
    nhood = dmat <= vdw_cutoff
    if npy.any(nhood):
        ats = npy.where(nhood)[0]
        for at in ats:
            parms = mol[at][0:3]
            probep = probe[0:3]

            Rmin = parms[1] + probep[1]
            E = npy.sqrt(parms[2]*probep[2])

            s = (Rmin/dmat[at])**6

            V += E*(s**2 - (2*s))

    return V

def calcElec(probe, mol, dmat):
    E = 0
    nhood = dmat <= elec_cutoff
    if npy.any(nhood):
        ats = npy.where(nhood)[0]
        for at in ats:
            q1 = mol[at][0]
            q2 = probe[0]

            E += (coul_k*q1*q2)/(dmat[at])
            
    return E


if __name__ == "__main__":
    import Grid as gr
    import sys
    
    if len(sys.argv) < 3:
        sys.exit("USAGE: python generateMIP.py PDB PROBE1 PROBE2 ...")

    pdbfile = sys.argv[1]

    # Load PDB
    pdb = bi.PDBModel(pdbfile)
    
    # Obtain parameters
    molecule = Molecule(pdb)
    
    # Build a grid over the molecule
    grid = gr.createAroundMolecule(pdb, 1., 4.)
    si,sj,sk = grid.shape
    size = grid.data.size
    
    # Choose a probe
    # Specified in the arguments??
    probes = sys.argv[2:]
    
    for probe in probes:
        
        params = probeparms[probe]
        grid.data *= 0  # Reset data to zeros
        
        n = 0
        for i in range(si):
            for j in range(sj):
                for k in range(sk):
                    n += 1
                    s = (float(n) / size)*100
                    if s in range(0,100,10):
                        print "%.2f %%\r"%(s)
                    xyz = grid.toCartesian((i,j,k))
                    probe_pos = params+xyz.tolist()
                    
                    dmat = distanceMat(molecule.xyz, xyz)
                    
                    vdw = calcVdW(probe_pos, molecule, dmat)
                    elec = calcElec(probe_pos, molecule, dmat)
                    
                    grid.data[i,j,k] = vdw + elec
                    
        grid.writeDX('%s_grid_test.dx'%probe)
        
    print "DONE"