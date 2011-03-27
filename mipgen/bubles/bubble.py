'''
Created on Mar 25, 2011

@author: psfriso
'''
import numpy as np
import Grid as gr
import sys
import BubFunctions as B

'''Initializating'''
print "  Bubles initilizating ..."

try:
    sys.argv[2]
    numout=B.check_cut(sys.argv[2])
except:
    numout=10

try:
    sys.argv[1]
    arg=sys.argv[1]
except:
    sys.exit("Grid file as first argument is mandatory")

file=B.check_file(arg)     

if file is None: sys.exit("Wrong Grid File")

fileout=open("/Users/psfriso/Desktop/Bubbles.out","w")

print "     Starting process @ ", B.time_stamp()

''' lectura de datos : uso Grid.py de dani'''
#infiere el formato
mygrid=gr.read(file)
origen=mygrid.origin
shape=mygrid.shape
size=mygrid.size
delta=mygrid.delta
#xplor = re.compile('.*xplor$',re.IGNORECASE)
#dx = re.compile('.*dx$',re.IGNORECASE)

'''numero de puntos candidatos a buenos
#dependera del tamanio de la molecula
# falta calibrar, cerca del 1%-0.5% del numero total
# parece ser suficiente'''
numtrial=15
''' numero de puntos con el que buscar vecinos, y sus energias
# cuanto mas grande mejor.'''
vecinos=10*numtrial
 
''' Cambiamos la forma a la grid'''
list=[]
for i in mygrid.data:
    list.append(i) 
''' intentamos prevenir irregularidades
#tambien se discretizan las energias
#las mallas no son precisas(poca densidad de puntos)
#es mas facil trabajar con valores suavizados'''
list=B.flat(B.limitar(list,size,shape),shape)
''' pasamos a formato {x,y,z,f(x,y,z)}'''
quedan=B.nonzero(list,shape)
hypsuper=B.toXYZ(origen,delta,shape,list,quedan)
'''ordenamos los valores por orden decreciente de energia'''
indices=hypsuper[:,3].argsort()

''' EXTRACT 20 MOST FAVORABLE POINTS'''
minima=np.zeros(shape=(numtrial,3))
for i in range(len(minima)):
    minima[i]=hypsuper[indices][i][0:3]
''' interpolando puntos. Haciendo mas denso el espacio importante'''
newmin=B.interpolation(minima,delta)
''' 
# List 100 most favorables
# doble lista con el objetivo de ajustar el radio de los puntos anteriores
# que son los puntos que inicialmente tienen valores minimo de energia
# los veinte primeros de ellos'''
minimaplus=np.zeros(shape=(vecinos,4))
for i in range(vecinos):
    minimaplus[i]=hypsuper[indices][i]

def core(newmin,minima,minimaplus,numout):
    '''# getting possible gauss'''
    pos_gauss=np.zeros(shape=(len(newmin),5))
    cont5=0
    for i in newmin:
        rad=B.max_coord(B.avg_5(B.expansion(i,minimaplus)))
        if (rad > 3.3):
            pass
        else:
            pos_gauss[cont5]=[rad,B.energy_inside(rad,i,minimaplus),i[0],i[1],i[2]]
            cont5+=1
    ''' ordeno las gausianas por la energia que abarcan dentro de 2*sd'''
    new_gauss=pos_gauss[np.argsort(pos_gauss[:,1])][0:int(numout)]
    return new_gauss

new_gauss=core(newmin,minima,minimaplus,numout)

''' Printing final results'''
print"        Working on file",file
print "           Searching for first",numout,"gaussians"
print "%14s%8s%9s%14s%6s" % ("X","Y","Z","ENERGY","SD")
print >> fileout, "Buble @",B.time_stamp()
print >>fileout,"%5s%8s%9s%14s%6s" % ("X","Y","Z","ENERGY","SD")
print >> fileout, " Gaussians generated for grid @",file
for i in new_gauss:
    print >>fileout,"%7.3f  %7.3f %7.3f  %10.2f  %6.3f" % (i[2],i[3],i[4],i[1],0.5*i[0])
    print "%14.3f  %7.3f %7.3f  %10.2f  %6.3f" % (i[2],i[3],i[4],i[1],0.5*i[0])
print '\n',"Output file Bubles.out was generated at",B.time_stamp()
print  "DONE :)"
