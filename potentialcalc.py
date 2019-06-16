# -*- coding: utf-8 -*-
"""
Created on my insomnia 

@author: Alarcon C.
"""
import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
#import time

DataName = input('Name of the input file(ej. "Data.dat"): ')
UnitsSpace = input('Units for the space(ej. "m", "cm", "mm", "km", "meter", "milimeter"): ')
UnitsPotential = input('Units for the potential(ej. "V", "mV", "kV", "Volt", "milivolt"): ')
NamePotData = input('Name for the output file(ej. "simulation.dat"):')

F = open(DataName)
argumentos = F.readline()
F.close()

#The first line of the imput file has the space configuration
#the first 2 values belong to the x axe, the next 2 belong to the y axe
#delta is the step used to discretize the axes, and tolerancia is the error accepted
X0, X1, Y0, Y1, delta, tolerancia = np.fromstring(argumentos, dtype=float, sep=' ')

Datos = np.loadtxt(DataName, skiprows = 1)

coord_xy=Datos[:,:2]
coord_xy[:,0]=coord_xy[:,0]-(X0)
coord_xy[:,1]=coord_xy[:,1]-(Y0)
coord_xy=coord_xy/delta
coord_xy[:, 1] = -coord_xy[:, 1] + abs(X1-X0)/delta
coord_xy.round().astype(int)
tmp = deepcopy(coord_xy)
coord_xy[:, 0] = tmp[:, 1]
coord_xy[:, 1] = tmp[:, 0]
coord_xy = tuple(map(tuple, coord_xy))

PotCte = np.zeros((np.arange(X0, X1+delta, delta).size, np.arange(Y0, Y1+delta, delta).size )) 
Mask = np.ones((np.arange(X0, X1+delta, delta).size, np.arange(Y0, Y1+delta, delta).size ))

for i in range(len(coord_xy)):
    PotCte[coord_xy[i]]=Datos[i, 2]
    Mask[coord_xy[i]] = 0.0

#plt.imshow(PotCte,cmap=plt.cm.jet)
#plt.colorbar()
#plt.show()

VPrevious = np.ones((np.arange(X0, X1+delta, delta).size, np.arange(Y0, Y1+delta, delta).size ));
VPrevious = VPrevious*Mask
VPrevious = VPrevious + PotCte

a, b  = VPrevious.shape
 
V1=deepcopy(VPrevious)
V2=deepcopy(VPrevious)
V3=deepcopy(VPrevious)
V4=deepcopy(VPrevious)

V1[:, 0:b-1]=V1[:,1:b]
V2[:, 1:b]  =V2[:,0:b-1]
V3[0:a-1,: ]=V3[1:a,:]
V4[1:a,: ]  =V4[0:a-1,:]

iteracion = 0
finish = 0

print("Start the calculation")

#t0 = time.clock()
while(finish == 0):
    iteracion = iteracion + 1
    VCurrent = (V1 + V2 + V3 + V4)*0.25
    VCurrent = VCurrent*Mask
    VCurrent = VCurrent + PotCte

    error = abs(VCurrent - VPrevious)
    #print('Error: ', error.max())
    V1=deepcopy(VCurrent)
    V2=deepcopy(VCurrent)
    V3=deepcopy(VCurrent)
    V4=deepcopy(VCurrent)
    
    V1[:, 0:b-1]=V1[:,1:b]
    V2[:, 1:b]  =V2[:,0:b-1]
    V3[0:a-1,: ]=V3[1:a,:]
    V4[1:a,: ]  =V4[0:a-1,:]

    VPrevious = deepcopy(VCurrent)
    if (tolerancia) > (error.max()):
        finish = 1

#t1 = time.clock()

np.savetxt(NamePotData,VCurrent,fmt='%.4e')#guardar imagen

del(VPrevious)
del(VCurrent)
del(V1)
del(V2)
del(V3)
del(V4)
del(PotCte)
del(Mask)
####################
####################

print (" Graphs \n\n")

MatPot = np.loadtxt(NamePotData)

PixelePerUnit= 1.0/delta
#m, n = MatPot.shape
length_y = Y1-Y0
length_x = X1-X0
#print("Image: ", length_y, "high" , UnitsSpace , length_x, UnitsSpace, "width")

x = np.linspace(X0, X1, np.arange(X0, X1+delta, delta).size);
y = np.linspace(Y1, Y0, np.arange(Y0, Y1+delta, delta).size);
X, Y = np.meshgrid(x,y); 

Ey, Ex = np.gradient(MatPot)
E = np.sqrt(Ex*Ex + Ey*Ey)

#Graph V
fig01=plt.figure()
plt.imshow(MatPot, extent=[X0, X1, Y0, Y1],cmap=plt.cm.jet)
plt.title('Potential')
plt.xlabel(UnitsSpace)
plt.ylabel(UnitsSpace)
cbar=plt.colorbar()
cbar.set_label(UnitsPotential)
plt.show(fig01)

#Graph E
fig01=plt.figure()
plt.imshow(E, extent=[X0, X1, Y0, Y1],cmap=plt.cm.jet)
plt.title('Magnitude Electric Field')
plt.xlabel(UnitsSpace)
plt.ylabel(UnitsSpace)
cbar=plt.colorbar()
cbar.set_label(UnitsPotential+'/'+UnitsSpace)
plt.show(fig01)

#Graph V cont 
fig01=plt.figure()
cs1 = plt.contourf(X, Y, MatPot, 25, cmap=plt.cm.jet)
cs2 = plt.contour(X, Y, MatPot, cs1.levels, colors='k')
plt.grid()
cbar=plt.colorbar(cs1)
cbar.set_label(UnitsPotential)
plt.xlabel(UnitsSpace)
plt.ylabel(UnitsSpace)
plt.title('Potential, contour lines')
plt.show(fig01)

#Graph V cont w/labels
fig01=plt.figure()
cs1 = plt.contourf(X, Y, MatPot, 25, cmap=plt.cm.jet)
cs2 = plt.contour(X, Y, MatPot, cs1.levels, colors='k')
plt.grid()
cbar=plt.colorbar(cs1)
cbar.set_label(UnitsPotential)
plt.clabel(cs2, colors='k', inline=1, fontsize=10)
plt.xlabel(UnitsSpace)
plt.ylabel(UnitsSpace)
plt.title('Potential, contour lines with labels')
plt.show(fig01)

#Graph E cont
fig01=plt.figure()
cs1 = plt.contourf(X, Y, E, 25, cmap=plt.cm.jet)
cs2 = plt.contour(X, Y, E, cs1.levels, colors='k')
plt.grid()
cbar=plt.colorbar(cs1)
cbar.set_label(UnitsPotential+'/'+UnitsSpace)
plt.xlabel(UnitsSpace)
plt.ylabel(UnitsSpace)
plt.title('Magnitude Elctric field, contour lines ')
plt.show(fig01)

#Graph E cont w/labels
fig01=plt.figure()
cs1 = plt.contourf(X, Y, E, 25, cmap=plt.cm.jet)
cs2 = plt.contour(X, Y, E, cs1.levels, colors='k')
plt.grid()
cbar=plt.colorbar(cs1)
cbar.set_label(UnitsPotential+'/'+UnitsSpace)
plt.clabel(cs2, inline=1, fontsize=9)
plt.xlabel(UnitsSpace)
plt.ylabel(UnitsSpace)
plt.title('Magnitude Elctric field, contour lines whith labels')
plt.show(fig01)

#Vector lines of E and V
fig01 = plt.figure()
cs01=plt.streamplot(X,Y,(-1)*Ex,Ey,
                    color = E,
                    cmap=plt.cm.jet,
                    linewidth=2,
                    arrowstyle='->',
                    arrowsize=1.5)
cs02 = plt.contour(X, Y, MatPot, 25, colors='b')
plt.xlim(X0, X1)
plt.ylim(Y0, Y1)
plt.xlabel(UnitsSpace)
plt.ylabel(UnitsSpace)
plt.title('Potential and  Vector electric field')
plt.show(fig01)



