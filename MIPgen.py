#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="dalvarez"
__date__ ="$04-feb-2011 18:39:53$"

import numpy as npy
import os
import string
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
    parser.add_option("-L", "--lib",dest="lib",default="parm/probes.lib",
            help="File containing the probes and its parameters. Default: parm/probes.lib")
#    parser.add_option("-A", "--parm",dest="parm",default=False,
#            help="File containing extra atom types missing in AmberFF or not identified automatically.")
    parser.add_option("-E","--eps",dest="eps",type="float", default=0.,
            help="Relative permitivity for the electrostatic calculations (float). If \
            0. is given, a Distance Dependent Relative Permitivity is used (default)")
    parser.add_option("-v","--vdw",dest="vdw",type="float",default=10.,
            help="VdW calculations Cutoff (float). Default: 10A")
    parser.add_option("-e","--elec",dest="elec",type="float",default=20.,
            help="Electrostatic calculations cutoff (float). Default: 20A")
    parser.add_option("-l","--list",default=False, action="store_true", dest="list",
            help="List available probes")

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
    "Function that calculates the van der Waals contribution"
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
    "Function to calculate the electrostatic contribution"
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

def getProbesFromFile(filename):
    """Reads a txt file like parm/probes.lib and returns a dictionary with the probes like:
    {NAME:[charge, vdwradii, vdweps, description],}"""
    import re

    # Read file
    # Format of the probes lines: NAME  CHARGE  VDWradii VDWeps  #DESCRIPTION
    f = open(filename,'r')
    lines = f.readlines()
    f.close()

    # Match lines beggining with # to skip them
    skip = re.compile(r'^#')

    result = {}
    
    for line in lines:
        if skip.match(line): continue
        else:
            line = line.split()
            result[line[0]] = map(float, line[1:4])
            result[line[0]].append(string.join(line[4:], ' '))

    return result

if __name__ == "__main__":
    from modules import Grid as gr
    from modules.amber import *
    import sys
    import os
    
    # WELCOME MESSAGE
    print '-'*90
    print '-'*20 + ' '*22 + 'MIPGEN' + ' '*22 + '-'*20
    print '-'*90
    print "MIPGEN is a python program that will calculate Molecular Interaction Potential grids"
    print "over a given molecule, that could be either a protein or a small organic compound (drug)."
    print
    print "The output will be a serie of grids with DX format (*.dx) that the user will be able"
    print "to visualize using any Molecular visualization program like VMD, PyMol, Chimera..."
    print
    print "For more information on dependencies and usage, please read the Documentation."
    print '-'*90
    print
    
    # Define options parser and get the options
    parser = defineParser()
    (options, args) = parser.parse_args(sys.argv[1:])

    # Load probes from lib
    probeparms = getProbesFromFile(options.lib)
    
    # Print probe list if options.l
    if options.list:
        print "Available probes:"
        for probe in probeparms.keys():
            print "\t+ "+probe+"\t\t"+probeparms[probe][-1]
        print
        print '-'*90
        sys.exit(0)
        
    if len(sys.argv[1:]) < 2:
        parser.error("Incorrect number of arguments. To get some help type -h or --help.")
        
    # Before beggining to work, check if AMBERTOOLS is installed
    if not os.environ.has_key('AMBERHOME'):
        sys.exit("ERROR: AMBERHOME environ path not found. Please check that you have a working installation of AMBERTOOLS.")
        
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

    print
    print "DONE"
    print
    print '-'*90