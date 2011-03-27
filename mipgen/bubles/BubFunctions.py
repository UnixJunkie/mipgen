'''
Created on Mar 26, 2011

@author: psfriso
'''
#def __init__(self, file=None, format=None, origin=None, shape=None,
#                spacing=None, size=None, value=0):
import numpy as np
import sys
import math    
import time
import os.path

def check_file(file):
    try:
        if(os.path.exists(file)):
            return file
    except IOError as e:
        print("({})".format(e))
        return None
#except IOError as e:
#    print("({})".format(e))

def check_cut(numout):
    try:
        numout=int(numout)
        if(numout):        
            numout=10
            print 'Incorrect numer of gaussians, adjusting to 10 (default)'
    except:
        numout=10
        print 'Incorrect numer of gaussians, adjusting to 10 (default)'
    return numout
def limitar(list,size,shape):
    cuantas_anula=0
    cuantas_hace=0
    shapei=shape[0]
    shapej=shape[1]
    shapek=shape[2]
    min_list=np.amin(list)
    lowercut=min_list*0.9
    uppercut=-0.2
    #print lowercut, uppercut, rango
    for i in range(shapei):
        for j in range(shapej):
            for k in range(shapek):
                cuantas_hace=cuantas_hace+1
                if (list[i][j][k]>uppercut or list[i][j][k]<lowercut):
                    list[i][j][k]=0
                    cuantas_anula+=1
    if(cuantas_hace!=size):
        sys.exit("Error @ masking list")
    return list



#fortiza las energias recibidas
def flat(list,shape):
    shapei=shape[0]
    shapej=shape[1]
    shapek=shape[2]
    uppercut=-0.2
    min_list=np.amin(list)
    lowercut=min_list*0.9
    rango=uppercut-lowercut
    interval=rango*0.2
    cont=0
# hago en lugar de un loop 5 veces el mismo bloque 
# se trata de un triple loop. esto puede ponerse muuuy lento con la estructura actual
# mejor usar los procesadores por separado...
 
    l=1
    lim1=uppercut-interval*(l-1)
    lim2=uppercut-interval*l
    for i in range(shapei):
        for j in range(shapej):
            for k in range(shapek):
                if (list[i][j][k]<= lim1 and list[i][j][k]>lim2):
                    list[i][j][k]=0.5*(lim1+lim2)
                    cont+=1

                    
    l=2
    lim1=uppercut-interval*(l-1)
    lim2=uppercut-interval*l
    for i in range(shapei):
        for j in range(shapej):
            for k in range(shapek):
                if (list[i][j][k]<= lim1 and list[i][j][k]>lim2):
                    list[i][j][k]=0.5*(lim1+lim2)
                    cont+=1

    l=3
    lim1=uppercut-interval*(l-1)
    lim2=uppercut-interval*l
    for i in range(shapei):
        for j in range(shapej):
            for k in range(shapek):
                if (list[i][j][k]<= lim1 and list[i][j][k]>lim2):
                    list[i][j][k]=0.5*(lim1+lim2)
                    cont+=1
                
    l=4
    lim1=uppercut-interval*(l-1)
    lim2=uppercut-interval*l
    for i in range(shapei):
        for j in range(shapej):
            for k in range(shapek):
                if (list[i][j][k]<=lim1 and list[i][j][k]>lim2):
                    list[i][j][k]=0.5*(lim1+lim2)
                    cont+=1

    l=5
    lim1=uppercut-interval*(l-1)
    lim2=uppercut-interval*l
    for i in range(shapei):
        for j in range(shapej):
            for k in range(shapek):
                if (list[i][j][k]<= lim1 and list[i][j][k]>lim2):
                    list[i][j][k]=0.5*(lim1+lim2)
                    cont+=1
    return list

# calcula la distancia euclideana de x y
def euclidean(x, y):
    if len(x) != len(y):
        raise ValueError, "vectors must be same length"
    suma=0
    for i in range(len(x)):
        suma += (x[i]-y[i])**2
    return math.sqrt(suma)

# pasa a {x,y,z,f(x,y,z)} los valores en dx
def toXYZ(origen,delta,shape,list,num):
    origenx=origen[0]
    origeny=origen[1]
    origenz=origen[2]
    dx=delta[0]
    dy=delta[1]
    dz=delta[2]
    shapei=shape[0]
    shapej=shape[1]
    shapek=shape[2]
    hypsuper=np.zeros(shape=(num,4))
    #energias=np.zeros(cuantas_quedan)
    cont2=0
    for i in range(shapei):
        for j in range(shapej):
            for k in range(shapek):
                if(list[i][j][k]<0):
                    hypsuper[cont2][0]=origenx+dx*i
                    hypsuper[cont2][1]=origeny+dy*j
                    hypsuper[cont2][2]=origenz+dz*k
                    hypsuper[cont2][3]=list[i][j][k]
#                    energias[cont2]=list[i][j][k]
                    cont2+=1

    if(num!=cont2):
        sys.exit("Error @ cambiar formato list")
    return hypsuper

#promedia los puntos en energias con sus vecinos
# es una forma de suavizar estas curvas horribles
def avg_5(vol_energia):
    long=len(vol_energia)
    r_avgener=np.zeros(shape=(long-2,2))
    for i in range(long-2):
        aux1=vol_energia[i][1]
        aux2=vol_energia[i+1][1]
        ncoord=vol_energia[i+1][0]
        aux3=vol_energia[i+2][1]
        avg=0.33333*(aux1+aux2+aux3)
        r_avgener[i]=[ncoord,avg]
    return r_avgener

# te devuelve la coordenada radial del maximo de energia englobada
def max_coord(radial):
    #print esfera
    long=len(radial)-1
    aux=radial[np.argsort(radial[:,1])][long][0]
    return aux

# calcula la energia de interaccion encerrada dentro de las 2 primeras
# desviaciones estandar de las gausianas
def energy_inside(radio,minimo,listaminimos):
    sumaener=0
    cont3=0
    long=len(listaminimos)
    for j in range(long):
        dist=euclidean(minimo,listaminimos[j][0:3])
        if(dist<=radio):
            sumaener+=listaminimos[j][3]
            cont3+=1
    cont3=0
    return sumaener

# vamos haciendo la gausiana cada vez mas grande para decidir donde cortar
def expansion(minimo_base,minimaplus):
    radioini=1.0
    radiofin=6.0
    radioincr=0.005*(radiofin-radioini)
    radio=radioini
    cont4=0
    vol_energia=np.zeros(shape=((radiofin-radioini)/radioincr+1,2))
    while(radio<radiofin):
        radio=radioini+radioincr*cont4
        energia=energy_inside(radio,minimo_base,minimaplus)
        norma=radio**3
        valor=energia**2/norma
        vol_energia[cont4]=[radio,valor]
        cont4+=1
        #print "div",valor,"ener",energia,"norma",norma,"r",radio
    #print vol_energia
    return vol_energia

# cuantos elementos no zero nos quedan. lo usaremos como control de no perder puntos
def nonzero(list,shape):
    shapei=shape[0]
    shapej=shape[1]
    shapek=shape[2]
    cont=0
    for i in range(shapei):
        for j in range(shapej):
            for k in range(shapek):
                if (list[i][j][k]!=0):
                    cont+=1            
    return cont

# interpola entre la malla para buscar posiciones candidatas a mejor centro de gaussiana
# es necesario para que podamos trabajar desde un buen principio con pocos minimos 
# se podria iterar varias veces...
def interpolation(lista,delta):
    entra=0
    long=len(lista)
    newlist=np.zeros(shape=(2*long,3))
    for i in range(long-1):
        for j in range(i+1,long):
            dist=euclidean(lista[i],lista[j])            
            if (dist <min(delta)*math.sqrt(2)):
                x=0.5*(lista[i][0]+lista[j][0])
                y=0.5*(lista[i][1]+lista[j][1])
                z=0.5*(lista[i][2]+lista[j][2])
                #print "entra",x,y,z,j,i
                newlist[entra]=[x,y,z]
                entra+=1
    for i in newlist:
        if(any(i)):
            lista=np.vstack([lista,i])
    return lista
def time_stamp():
# Using datetime.strptime()
    now = time.localtime(time.time())
    return time.strftime("%d/%m/%y %H:%M", now)