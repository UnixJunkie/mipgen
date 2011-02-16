# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="dalvarez"
__date__ ="$04-feb-2011 18:40:50$"

import re
import os, sys
import os.path as path
import time
import fnmatch

import numpy as npy

try:
    import Biskit as bi
    _biskit = True
except:
    _biskit = False

try:
    import Bio.PDB as bpdb
    _bio = True
except:
    _bio = False

if not _bio and not _biskit:
    error("Biskit or BioPython need to be installed to use some of the \
    functionalities of the module. createAroundMolecule will not work")

_dxtemplate="""#Grid generated with GRID.py
# Time: %s
# Source: %s
object 1 class gridpositions counts %i %i %i
origin %.2f %.2f %.2f
delta %.5f   0   0
delta 0   %.5f   0
delta 0   0   %.5f
object 2 class gridconnections counts %i %i %i
object 3 class array type double rank 0 items %i data follows
"""

_xplortemplate="""XPLOR Format\n       1\nFree binding energies GRID generated with count2DG.py
%8d%8d%8d%8d%8d%8d%8d%8d%8d\n%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f
ZYX
"""

class Grid:
    
    def __init__(self, file=None, format=None, origin=None, shape=None,
                spacing=None, size=None, value=0):

        """
        Constructor method for a 3D grid

        Two ways of constructiong the instance:

        1 - Reading from a file with DX or XPLOR format

        If file is given, it will try to read from that file.
        Formats supported are DX and XPLOR. Format argument is mandatory

        2 - From scratch

        If file is None, the constructor will use the rest of parameters
        to build a 'value' filled grid instance with the specified origin and spacing.

        The size of the grids can be given in two formats aswell:

            2.1 - Using a shape: It is the number of grid points in each axis.
                                    Should be a tuple or list.
            2.2 - Using a size: It is the size in Angstroms in each axis.
                                    Witht the spacing information, the number of
                                    grid points will be then calculated.

        TYPES:

        file    -   String.
        format  -   String. At the moment only 'DX' or 'XPLOR' is supported.

        origin  -   List [x,y,z]. Origin in cartesian coordinates.
        spacing -   Float or List. This is the size of the cells in Angstroms.
                    If Float, all axis will have the same spacing.
                    If List, each axis can have different spacing parameter.
        shape   -   List with number of points in x,y,z axis.
        size    -   List or float.
                    List of length 3 with size of each side in Angstroms.
                    Float for same size all sides. Also in Angstroms.
        value   -   Int or Float. Filling value for all points in the grid.

        """
        # Constructor method 1
        if file and format:
            self.source = path.abspath(file)
            self._readfile(format)
        
        # Constructor method 2
        elif shape is not None or size is not None and not file:
            if isinstance(origin, list) or isinstance(origin, tuple): self.origin = map(float, origin)
            else: print "ERROR. Wrong format for origin parameter"

            if isinstance(spacing, float) or isinstance(spacing, int): self.delta = map(float, [spacing, spacing, spacing])
            elif isinstance(spacing, list): self.delta = map(float, spacing)
            else: print "ERROR. Wrong format for spacing parameter"

            # Size or shape to
            # create numpy array in data
            if size is not None:
                shape = npy.round(npy.array(size) / self.delta)

            self.data = npy.zeros(shape)
            self.data += value

            self.source = "User created Grid"

        # Not enough arguments
        else:
            error("Could not create the grid instance. Possible reasons: File or format missing. Neither size nor shape specified.")

        # Assign some useful attributes:
        self.shape = self.data.shape
        self.size = self.data.size
        
    
    def _readfile(self, format):
        "Use _readDX or _readXPLOR depending on the format"
        
        if format == 'DX':
            self._readDX()
        elif format == 'XPLOR':
            self._readXPLOR()
        else:
            error("File format not valid. Only DX and XPLOR formats supported.")


    def _readDX(self):
        "Read DX format grid files"

        f=open(self.source,"r")
        
        #read the header
        header=""
        for i in range(3) : header= header + f.readline()
        
        #read the grid size
        r=re.compile('\w+')
        gsize=r.findall(f.readline())
        gsize=[int(gsize[-3]),int(gsize[-2]),int(gsize[-1])]

        #read the origin of the system
        line=f.readline().split()
        origin=[float(line[-3]),float(line[-2]),float(line[-1])]

        #read grid spacing
        line=f.readline().split()
        deltax=[float(line[-3]),float(line[-2]),float(line[-1])]
        line=f.readline().split()
        deltay=[float(line[-3]),float(line[-2]),float(line[-1])]
        line=f.readline().split()
        deltaz=[float(line[-3]),float(line[-2]),float(line[-1])]

        #pay attention here, this assumes always orthogonal normalized space, but normally it should be ok
        delta=[deltax[0],deltay[1],deltaz[2]]
        
        #read the number of data
        f.readline()
        r=re.compile('\d+')
        n_entries=int(r.findall(f.readline())[2])

        #check correpondence with expected data points
        if(n_entries!=gsize[0]*gsize[1]*gsize[2]) :
            sys.exit("Error reading the file. The number of expected data points \
            does not correspond to the number of labeled data points in the header.")

        #load data into numpy array
        #reshaping to fit grid format (it keeps Z fast, Y medium, X slow data organization)
        grid = npy.fromstring(f.read(), sep=' ', dtype=float).reshape(gsize)

        if grid.size != n_entries:
            sys.exit("Error reading the file. The number of expected data points\
            does not correspond to the number of labeled data points in the header.")
        f.close()
        
        self.data= grid
        self.origin = origin
        self.delta = delta

    def _readXPLOR(self):
        "Read XPLOR format grid files"

        #Open XPLOR Grid File
        f = open(self.source, 'r')

        #Skip 1st line
        f.readline()
        
        #Header lines:
        head = f.readline().split()
        head = int(head[0])
        while head:
            f.readline()
            head-=1
            
        #Read grid format
        format1 = f.readline().split()
        for i in xrange(len(format1)): format1[i] = int(format1[i])
        nx, minx, maxx, ny, miny, maxy, nz, minz, maxz = format1

        format2 = f.readline().split()
        for i in xrange(len(format2)): format2[i] = float(format2[i])
        dx, dy, dz, deltax, deltay, deltaz = format2

        #Skip ZYX line
        f.readline()

        #Calculate size to read
        #XPLOR Format has ZYX format, X fast, Y medium, Z slow.
        #X data is sorted in 6 columns
        #Each data is 12 bits and each line is 1 bit more (newline charater)
        #This format is repeated ny times. Then a new Z block starts (one line to be skipped then).
        size = int((nx*12)+npy.ceil(nx/6.))
        grid = npy.zeros([nx,ny,nz])
        
        for z in xrange(nz):
            f.readline() #skip block identifier line
            for y in xrange(ny):
                grid[:,y,z]=readBlock(f, size, 'float')
                
        #Prepare standard descriptors:
        delta = [(dx/nx), (dy/ny), (dz/nz)]
        origin = [minx*delta[0],miny*delta[1],minz*delta[2]]
        
        #Assign attributes
        self.data = grid
        self.origin = origin
        self.delta = delta

    def toCartesian(self, index):
        """Conversion from grid index to cartesian coordinates.

        index must be a list or tuple of length 3.
        
        If index is not valid, None will be returned.
        Numpy array will be returned if all is correct.
        """

        if len(index) != 3:
            return None
        
        point = npy.array(index)
        valid = npy.all((point < self.shape) * (point >= 0))

        if not valid:
            return None
        
        return self.origin + point * self.delta
                
        
    def toIndex(self, coords):
        """Conversion from coordinates to grid index if fits inside.
        Return None otherwise.Coords has to be tuple or list with length 3.
        Numpy array will be returned if all is correct.
        """

        if len(coords) != 3:
            return None

        point = npy.array(coords)
        maxcoord = npy.array(self.shape) * self.delta + self.origin
        valid = npy.all((point >= self.origin) * (point < maxcoord))
        
        if not valid:
            return None

        return npy.floor((point - self.origin) / self.delta).astype(int)
        
        
    def writeXPLOR(self, xplorname):
        """
        Write data into XPLOR formated file.
        """
        grid = self.data
        origin = self.origin
        delta = self.delta

        nx, ny, nz = grid.shape

        mind = npy.round(npy.array(origin) / delta)
        maxd = npy.round(mind + grid.shape - 1)
        space = npy.array(delta)*grid.shape         # Calculate total system expansion (distance) in x, y, z.
        angle = 90.                                 # Assuming rectangular shape

        of = open(xplorname, 'w')
        #Header and grid format
        of.write(_xplortemplate%(nx, mind[0], maxd[0], ny, mind[1], maxd[1], nz, mind[2],
        maxd[2], space[0], space[1], space[2], angle, angle, angle))

        #Write Grid data
        for z in xrange(nz):
            of.write("%8d\n"%(mind[2]))
            for y in xrange(ny):
                nl = 0
                for el in grid[:,y,z]:
                    nl+=1
                    of.write("%12.5f"%(el),)
                    if nl == 6: of.write("\n"); nl=0
                if (nx%6 != 0): of.write("\n")
            mind[2]+=1
        of.close()
            


    def writeDX(self, dxname):
        """
        Writes data into a DX Formated File
        """
        grid = self.data
        origin = self.origin
        delta = self.delta
        dxf = open(dxname,'w')
        dxf.write(_dxtemplate%(time.ctime(),self.source,grid.shape[0],
                                            grid.shape[1],grid.shape[2],
                                            origin[0],origin[1],origin[2],
                                            delta[0],delta[1],delta[2],
                                            grid.shape[0],grid.shape[1],
                                            grid.shape[2],grid.size))

        count=0
        for data in grid.flat:
            if count==3: dxf.write('\n'); count=0
            dxf.write('%6.3f\t'%(float(data)))
            count+=1
        dxf.write('\n')
        dxf.close()


def error(message):
    print >> sys.stderr, "ERROR: ",message

def readBlock(fh, size, dType):
    s = npy.fromstring(fh.read(size),sep=" ",dtype=dType)
    return s.copy()

def read(file, format=None):
    "Helper function for constructing a Grid instance from a valid file. \
    It reads DX and XPLOR format grids."

    xplor = re.compile('.*xplor$',re.IGNORECASE)
    dx = re.compile('.*dx$',re.IGNORECASE)
    
    # Check wheter file exists or not
    if path.exists(file):

        # Check wheter format is given or not
        # Then try to match format names
        if not format:
            if dx.match(file):
                format = 'DX'
            elif xplor.match(file):
                format = 'XPLOR'
            else:
                error("Format not recognized. Provide the \
                        format of the file as argument format='DX' or change \
                        the filename extension to DX or XPLOR for autodetection.")

                return None
        else:
            if not dx.match(format) and not xplor.match(format):
                # Format unrecognized
                error("Format not valid. Should be 'dx' or 'xplor'.")
                return None

        return Grid(file=file, format=format)

    else:
        # File does not exist
        error("File path not correct or file does not exist.")
        return None

def create(**kwargs):
    """Helper function which redirects to the creation of the instance.
    Just to keep nice names"""
    return Grid(**kwargs)


def createAroundMolecule(molecule, spacing, buffer):
    """Build a Grid instance around molecule with the given spacing and
    some buffer around.
    
    Molecule is a string of a PDB file or
    an instance of a Structure in Biopython or
    a PDBModel instance in Biskit.

    Molecule can even be a numpy array with the coordinates.

    buffer is Float or int
    """
    
    if isinstance(molecule, npy.ndarray):
        if molecule.shape[1] == 3:
            xyz = molecule
    else:
        xyz = coordsFromMolec(molecule)

    if xyz is None:
        error("Error obtaining coordinates for the molecule")

    # We have the coordinates
    # let's see what is the minimum and maximum coordinates
    # stablish the distance and origin according to the buffer
    min = xyz.min(axis=0)
    max = xyz.max(axis=0)
    size = (max - min) + (buffer*2)
    origin = min - buffer
    
    return Grid(size=size, origin=origin.tolist(), spacing=spacing)

def coordsFromMolec(molecule):
    """Molecule is a string of a PDB file or
    an instance of a Structure in Biopython or
    a PDBModel instance in Biskit.

    This function loads the file and/or gets the xyz coords
    from the instances created"""
    xyz = None
    
    if _bio:
        if isinstance(molecule, bpdb.Structure.Structure):
            xyz = getStructureCoords(molecule)

    if _biskit and xyz is None:
        if isinstance(molecule, bi.PDBModel):
            xyz = molecule.xyz

    if isinstance(molecule, str) and xyz is None:
        if _biskit:
            xyz = bi.PDBModel(molecule).xyz
        elif _bio:
            parser = bpdb.PDBParser()
            molecule = parser.get_structure('mol',molecule)
            xyz = getStructureCoords(molecule)

    return xyz


def getStructureCoords(Struct):

    atoms = Struct.get_atoms()
    coords = [at.get_coord().tolist() for at in atoms]

    return npy.array(coords)


if __name__ == "__main__":
    pdb = bi.PDBModel('../cheachey.pdb')
    a = createAroundMolecule(pdb.xyz,0.5,4)
    a.writeXPLOR('test.xplor')
