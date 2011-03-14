#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="dalvarez"
__date__ ="$04-feb-2011 18:39:53$"

import sqlite3
import numpy as npy
from os.path import split, splitext, exists
import subprocess as sub
from optparse import OptionParser

def defineParser():
    " Define the Option parser"
    parser = OptionParser()
    
    parser.add_option("-p","--prot",dest="prot", default=False,
            help="Protein file (PDB format)")
    parser.add_option("-m","--molec",dest="molec", default=False,
            help="Molecule file (PDB format)")
    parser.add_option("-r", "--probe",dest="probes",action="append",
            help="Append probe names to calculate the MIPs. This flag can be used\
             more than once.")
    parser.add_option("-o","--out",dest="out",default="MIP",
            help="Output name prefix. Use some name descriptive of your job. Default: MIP")
    parser.add_option("-E","--eps",dest="eps",type="float", default=0.,
            help="Relative permitivity for the electrostatic calculations (float). If \
            0. is given, a Distance Dependent Relative Permitivity is used (default)")
    parser.add_option("-v","--vdw",dest="vdw",type="float",default=10.,
            help="VdW calculations Cutoff (float). Default: 10A")
    parser.add_option("-e","--elec",dest="elec",type="float",default=20.,
            help="Electrostatic calculations cutoff (float). Default: 20A")

    return parser

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

def calcElec(probe, mol, dmat, diel_const):
    econst = bool(diel_const)
    E = 0
    if econst: coul_k = 332./(diel_const)
    nhood = dmat <= elec_cutoff
    if npy.any(nhood):
        ats = npy.where(nhood)[0]
        for at in ats:
            q1 = mol[at][0]
            q2 = probe[0]
            r = dmat[at]
            if not econst: coul_k = 332./(4*r)
            E += (coul_k*q1*q2)/(r)

    return E

class Leaper:
    """Handle Leap I/O"""
    def __init__(self,  cwd=None):
        self.tleap = sub.Popen("tleap -f - ",  shell=True,   stdin = sub.PIPE,  stdout=sub.PIPE,  stderr=sub.PIPE, cwd=cwd)
        self._in = self.tleap.stdin
        self._out = self.tleap.stdout
        self._err = self.tleap.stderr
        
    def command(self,  command):
        """
        Command must be a string corresponding to a correct Leap command.
        Returns the output of leap for the given command.
        tLeap is a bit tricky to be controlled by stdin because it terminates when it reaches EOF when reading STDOUT.
        We use a fake command to stop reading STDOUT right before the EOF is reached.
        """
        self._in.write(command+'\n')
        self._in.write('lastCommandOut\n')
        out = []
        while 1:
            line = self._out.readline().strip()
            if 'ERROR: syntax error' in line: break
            if line: out.append(line)
        return out
        
    def close(self):
        self.tleap.terminate()
        del self.tleap  

class AmberMolecule:
    
    def __init__(self, pdbFile, mode):
        """Build an instance for the molecule containing
        for each atom, the coordinates and non-bonded parameters
        (mass, charge, vdw radii and epsilon)

        Input:  atom_names - string with names as in PDB
                xyz - numpy.ndarray with the coordinates
        """
        
        self.atoms = self.__getParamsXYZ(pdbFile, mode)
#        print self.atoms
#        sys.exit()
        self.natoms = len(self.atoms)
        self.xyz = npy.array([at[3:] for at in self.atoms])
        
        
    def __getParamsXYZ(self, pdbFile, mode):
        """
        Obtain for each atom the correspondent parameters for
        Vander Waals and Coulomb calculations along with their coordinates.
        This is done using an external software: AmberTools v.1.1.

        For organic molecules, GAFF force field should be used, and PARM99
        is used for proteins atom type identification.

        Thus two modes can be chosen in the mode parameter of the function:

        Mode = 'gaff'   --> Organic molecule. Treat with antechamber.
        Mode = 'parm'   --> Protein. Treat with tLeap and then assign Amber FF params.

        The function returns a list with this format:
        [[charge, vdw, eps, x, y, z],[..],..]
        """
        
        # If its a protein, 'parm' mode, first give a try to normalize a bit
        # the file using tLeap from AmberTools. It happens that the PDB can
        # have different occupancies and/or extrange names. This will filter a bit
        # this abnormalities
        if mode == 'parm':
            pdb = self.__fixPDBtLeap(pdbFile)
            if pdb: pdbFile = pdb
            else: pass # Will continue without pre-process and good luck!
        
        # GAFF is used for organic molecules. Antechamber is used to prepare them.
        # PARM99 is used for proteins. Parameters are obtained from amber.db
        if mode == 'gaff':
            # First obtain atom types and charges
            types_charges = self.__gaffParams(pdbFile)
            res_at_xyz = self.__getAtsXYZ(pdbFile)
            
            # Atoms in types_charges don't have the same order as in the PDB file
            # we have to match atom indices with the coordinates
            at_xyz_dict = dict(zip([at[1] for at in res_at_xyz],[at[2:] for at in res_at_xyz]))
            xyz = [at_xyz_dict[at[0]] for at in types_charges]
            
            # Assign corresponding VdW parameters using gaff.db
            params = self.__addVdWtoGaff(types_charges)
            
            return [params[i]+xyz[i] for i in range(len(params))]
        
        elif mode == 'parm':
            # First extract from the PDB the coordinates, atom names and residue
            # name for each atom
            res_at_xyz = self.__getAtsXYZ(pdbFile)
            res_at = [[at[0],at[1]] for at in res_at_xyz]
            params, miss = self.__amberParams(res_at)
            xyz = [at[2:] for at in res_at_xyz if at[0:2] not in miss]
            return [params[i]+xyz[i] for i in range(len(params))]

    def __queryDB(self, sqliteConnection, DB, attype=None, atname=None, resname=None):
        "Query database amber.db or gaff.db by atom name or by atom type and return result.\
        DB should be either 'amber' or 'gaff'."

        if DB == 'gaff':
            table = 'gafftypes'
        elif DB == 'amber':
            table = 'ambertypes'
        else:
            return False
        
        # Create a pointer
        c = sqliteConnection.cursor()

        # Prepare where clause:
        if atname and resname:
            where = "resname='%s' and atname='%s'"%(resname, atname)
        elif attype:
            where = "attype='%s'"%attype
            
        # Query it!
        c.execute("SELECT * FROM %s WHERE %s"%(table,where))
        row=c.fetchall()
        if len(row) != 1:
            # No results
            print 'missing atom'
            print "At:",atname,resname,attype
            return 'miss',resname,atname
        else:
            # Correct fetching
            return row[0]
            
    def __addVdWtoGaff(self,types_charges):
        "Add the corresponding VdW parameters from gaff.db to the given atomtypes"

        # Open sqlite connection and query the db for each atom type
        conn = sqlite3.connect('ffparm/gaff.db')
        # TODO if some atom is missing, catch the error!
        vdw=[self.__queryDB(conn, 'gaff', attype=at[1])[2:] for at in types_charges]
        conn.close()
        
        # Return a list Charge,VdWRadii,Eps
        return [[float(types_charges[i][2])]+list(vdw[i]) for i in range(len(types_charges))]

        
    def __fixPDBtLeap(self, pdbFile):
        "Pre-process PDB using AmberTools software. Just load it and save it again.\
        It will change some atom names and some residues hopefully fixing some strange things.\
        Will also add missing Hydrognes."
        outname = pdbFile.replace('.','_tleap.')
        tleap = Leaper()
        tleap.command('p = loadPdb %s'%pdbFile)
        tleap.command('savePdb p %s'%outname)
        tleap.close()
        
        if exists(outname):
            return outname
        else:
            return False
        
    def __getAtsXYZ(self,pdbFile):
        "Given a pdb file, return the list of residue names and atoms names for each atom."
        from PDBParser import readResAtCoordFromPDB
        molec = readResAtCoordFromPDB(pdbFile)
        return molec


    def __amberParams(self, res_at):
        """Returns for each pair residue_atomname a tuple (charge, radii, epsilon)
        Using Amber 99 FF"""
        # Open sqlite connection and query the db for each atom type
        conn = sqlite3.connect('ffparm/amber.db')
        # TODO if some atom is missing, catch the error!
        params = []
        miss = []
        for at in res_at:
            if at[1] == 'OXT': at[0] = 'C'+at[0]
            if at[0] == 'HIS': at[0] = 'HIE'
            p = self.__queryDB(conn, 'amber', resname=at[0], atname=at[1])
            if p[0] != 'miss': params.append(list(p)[4:])
            else: miss.append(p[1:])
        conn.close()
#        print params
        # Return a list Charge,VdWRadii,Eps
        return params, miss
        
    def __gaffParams(self, pdbFile):
        """Returns for each atom in pdbInstance a tuple (charge, radii, epsilon, x, y, z)
        Using Amber GAFF FF."""
        from amberprep import readAmberPrep
        
        pdbin = pdbFile
        prepin = splitext(split(pdbFile)[1])[0] + '.prep'
        ante_cmd = "antechamber -i %s -fi pdb -o %s -fo prepi -c bcc -pf"%(pdbin,prepin)
        
        print "Running antechamber to calculate charges and assign atomtypes..."
        if not exists(prepin):
            antechamber = sub.Popen(ante_cmd, shell=True, stdout=sub.PIPE, stderr=sub.PIPE)
            antechamber.wait()
#            antechamber.terminate()
        
        if exists(prepin):
            molec = readAmberPrep(prepin)
            # Will return only first residue as it should only contain 1 :)
            return molec[molec.keys()[0]]
        else:
            print antechamber.stdout.read()
            print antechamber.stderr.read()
            return None

    
    def __getitem__(self, idx):
        return self.atoms[idx]


if __name__ == "__main__":
    import Grid as gr
    import sys
    import os
    
    # PROBE Parameters (CHARGE, VDW RADII, EPS)
    probeparms = {'C':[0.0,1.95, 0.07], # CD from leucine
             'HDON':[0.446255,0., 0.0], # from hydrogen in SER (HO)
             'pos':[1,0, 0.0], # pure +1 charge
             'N+':[0.633325,1.85, 0.18], # taking N3 charge + 3H charges and N3 radii + a bit increment to account for Hs
             'neg':[-1,0, 0.0]} # pure -1 charge
             
    # Define options parser and get the options
    parser = defineParser()
    (options, args) = parser.parse_args(sys.argv[1:])
    
    if len(sys.argv[1:]) < 2:
        parser.error("Incorrect number of arguments. To get some help type -h or --help.")

    # Output name
    outprefix = options.out
    
    # Parse VdW and Elec options
    vdw_cutoff = options.vdw
    elec_cutoff = options.elec
    diel_const = options.eps
#    print "Using diel_const: ",diel_const, bool(diel_const)
    
    # Check if the user has given a molecule or a protein argument
    # They will be treated differently when assigning atom types
    if options.molec and options.prot:
        parser.error("Protein and Molecule options are exclusive. To get some help type -h or --help.")
    elif options.molec:
        pdbfile = os.path.exists(options.molec) and options.molec or False
        mode = "gaff"
    elif options.prot:
        pdbfile = os.path.exists(options.prot) and options.prot or False
        mode = "parm"
    else:
        parser.error("Either Protein or Molecule have to be given. To get some help type -h or --help.")

    if not pdbfile: parser.error("Path to PDB does not exists.")
    
    # Check that all the probes chosen are defined
    probes = options.probes
    [parser.error("Error in probe name %s. Check manual for available probes."%p)
        for p in probes if p not in probeparms.keys()]
    
    # Obtain parameters
    print "Obtaining parameters..."
    molecule = AmberMolecule(pdbfile,mode)
    
    # Build a grid over the molecule
    print "Building grid over the molecule..."
    grid = gr.createAroundMolecule(molecule.xyz, 1., 4.)
    si,sj,sk = grid.shape
    size = grid.data.size
    
    for probe in probes:
        print "\t+ Running probe: ",probe
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
                    elec = calcElec(probe_pos, molecule, dmat, diel_const)
                    
                    grid.data[i,j,k] = vdw + elec
                    
        grid.writeDX(outprefix+'_%s.dx'%probe)
        
    print "DONE"