import parser
#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="dalvarez"
__date__ ="$04-feb-2011 18:39:53$"

import Biskit as bi
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
    parser.add_option("-E","--eps",dest="eps",type="float", default=0.,
            help="Relative permitivity for the electrostatic calculations (float). If \
            0. is given, a Distance Dependent Relative Permitivity is used (default)")
    parser.add_option("-v","--vdw",dest="vdw",type="float",default=10.,
            help="VdW calculations Cutoff (float). Default: 10A")
    parser.add_option("-e","--elec",dest="elec",type="float",default=20.,
            help="Electrostatic calculations cutoff (float). Default: 20A")
    parser.add_option("-r", "--probe",dest="probes",action="append",
            help="Append probe names to calculate the MIPs. This flag may be used\
             once per probe.")
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

def calcElec(probe, mol, dmat, dist_dependent=False):
    E = 0
    if not dist_dependent: coul_k = 332./(dielectric_constant)
    nhood = dmat <= elec_cutoff
    if npy.any(nhood):
        ats = npy.where(nhood)[0]
        for at in ats:
            q1 = mol[at][0]
            q2 = probe[0]
            r = dmat[at]
            if dist_dependent: coul_k = 332./(4*r)
            E += (coul_k*q1*q2)/(r)

    return E

class Leaper:
    """Handle Leap I/O"""
    def __init__(self,  cwd=None):
        self.tleap = subprocess.Popen("tleap -f - ",  shell=True,   stdin = subprocess.PIPE,  stdout=subprocess.PIPE,  stderr=subprocess.PIPE, cwd=cwd)
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
        self.command("quit")
        self.tleap.terminate()
        del self.tleap  

class Molecule:

    def __init__(self, pdbFile, mode):
        """Build an instance for the molecule containing
        for each atom, the coordinates and non-bonded parameters
        (mass, charge, vdw radii and epsilon)

        Input:  atom_names - string with names as in PDB
                xyz - numpy.ndarray with the coordinates
        """
        
        self.atoms, self.missing = self.__getParamsXYZ(pdbFile, mode)
        self.natoms = len(self.atoms)
        self.xyz = npy.array([at[0][3:] for at in self.atoms])
        
        
    def __getParamsXYZ(self, pdbFile, mode):
        """
        Obtain for each atom the correspondent parameters for
        Vander Waals and Coulomb calculations along with their coordinates.
        This is done using an external software: AmberTools v.1.1.

        For organic molecules, GAFF force field should be used, and PARM99
        is used for proteins atom type identification.

        Thus two modes can be chosen in the mode parameter of the function:

        Mode = 'gaff'   --> Organic molecule. Treat with antechamber.
        Mode = 'parm'   --> Protein. Treat with Amber FF in tleap.

        The function returns a list with this format:
        [[charge, vdw, eps, x, y, z],[..],..]
        """
        # Prepare Prepin and read it with tleap or antechamber
        # depending on the input option.
        #
        # Prepin file contains the assigned atom types and charges
        # later, on the amber.db we will fetch the lacking VdW parameters
        #
        # GAFF is used for organic molecules. Antechamber is used to prepare them.
        # PARM99 is used for proteins. tLeap is used for preparing the prepin.
        if mode == 'gaff':
            res_types_charges = self.__gaffParams(pdbFile)
        elif mode == 'parm':
            res_types_charges = self.__amberParams(pdbFile)

        # The output of the previous step is a dictionary containing all the residues
        # identified in the prepin file and the atoms/charges per residue.
        #
        # Next step is to fetch the VdW parameters for the atomtypes assigned in the DB.
        #
        # We have to take also the coordinates for the atoms. So firstly will use Biskit module to read
        # the pdb file and get the coordinates.
        
        xyz = bi.PDBModel(pdbFile).xyz
        
        if res_types_charges:
            for res in res_types_charges.keys():
                res = res_types_charges[res]
                for at in range(len(res)):
                    pass

    def __amberParams(self, pdbFile):
        """Returns for each atom in pdbInstance a tuple (charge, radii, epsilon, x, y, z)
        Using Amber GAFF FF."""
        pdbin = pdbFile
        prepin = splitext(split(pdbFile)[1])[0] + '.prep'
        
        # Open a tleap connection
        # and send commands to build a prepin file from the pdb
        leap = Leaper()
        leap.command("source ffparm/parm99.dat")
        leap.command("p = loadPdb %s"%pdbFile)
        leap.command("saveAmberPrep p %s"%prepin)
        leap.close()
        
        # Check that prepin file was created
        # read it and return output
        if exists(prepin):
            return self.__readAmberPrep(prepin)
        else:
            return None
        
    def __gaffParams(self, pdbFile):
        """Returns for each atom in pdbInstance a tuple (charge, radii, epsilon, x, y, z)
        Using Amber GAFF FF."""
        pdbin = pdbFile
        prepin = splitext(split(pdbFile)[1])[0] + '.prep'
        ante_cmd = "antechamber -i %s -fi pdb -o %s -fo prepi -c bcc -pf"%(pdbin,prepin)

        print "Running antechamber to calculate charges and assign atomtypes..."
        antechamber = sub.Popen(ante_cmd, shell=True, stdout=sub.PIPE, stderr=sub.PIPE)
        antechamber.wait()

        if exists(prepin):
            return self.__readAmberPrep(prepin)
        else:
            return None

    def __readAmberPrep(self, prepfile):
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
    
    def __charmmParams(self,pdbFile):
        """Returns for each atom in pdbInstance a tuple (charge, radii, epsilon, x, y, z)
        Using charmm DB."""
        
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

        # Create Biskit.PDBModel Instance and
        # store coordinates, residue names and atom names
        # for faster access during the loop
        pdbInstance = bi.PDBModel(pdbFile)
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


if __name__ == "__main__":
    import Grid as gr
    import sys
    import os

    # PROBE Parameters (CHARGE, VDW RADII, EPS)
    probeparms = {'C':[0.0,1.95,-0.07], # CD from leucine
             'O':[-0.01,1.65,-0.12],
             'H+':[1,1,-0.1],
             'N':[-0.47,1.85,-0.2],
             'S':[-0.09,2.1,-0.47],
             'neg':[-1,2.5,-0.7]}

    # Define options parser and get the options
    parser = defineParser()
    (options, args) = parser.parse_args(sys.argv[1:])
    
    if len(sys.argv[1:]) < 2:
        parser.error("Incorrect number of arguments. To get some help type -h or --help.")
        
    # Parse VdW and Elec options
    vdw_cutoff = options.vdw
    elec_cutoff = options.elec
    diel_const = options.eps
    
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
    probes = options.parser
    [parser.error("Error in probe name %s. Check manual for available probes."%p)
        for p in probes if p not in probeparms.keys()]
    
    # Obtain parameters
    molecule = Molecule(pdbfile,mode)
    
    # Build a grid over the molecule
    grid = gr.createAroundMolecule(pdb, 1., 4.)
    si,sj,sk = grid.shape
    size = grid.data.size
    
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
                    elec = calcElec(probe_pos, molecule, dmat, dist_dependent=bool(diel_const))
                    
                    grid.data[i,j,k] = vdw + elec
                    
        grid.writeDX('%s_grid_test.dx'%probe)
        
    print "DONE"