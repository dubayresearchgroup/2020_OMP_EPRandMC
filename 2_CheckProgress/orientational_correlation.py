#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 14:57:25 2019

@author: katerilab
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 13:55:53 2019

@author: nhunguyen
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 09:53:07 2019

@author: katerilab
"""


import numpy as np 
import matplotlib.pyplot as plt
from readTraj_2D import ReadSnaps


def Periodic(distance, boxlength): # To calculate the distance with periodic boundary condition
    
    L = distance
    B = boxlength
    L = (L + B /2)%B - B/2
    return L

def Dist(p1, p2, boxL):

    dist = abs(p1 - p2)
    if dist > 0.5 * boxL:
        dist = boxL - dist

    return dist


def Type_DataMix(mol, typ, x_raw, y_raw, dimL, noWant):
    
    molID = mol[np.where(typ != noWant)]
    modTyp = typ[np.where(typ != noWant)]
    x = x_raw[np.where(typ != noWant)] 
    y = y_raw[np.where(typ != noWant)]

    return molID, modTyp, x, y

def Type_DataCenter(typ, x_raw, y_raw, dimL, want):
    
    modType = typ[np.where(typ == want)]
    x = x_raw[np.where(typ == want)] 
    y = y_raw[np.where(typ == want)]
    
    return x, y

def Type_Data(typ, x_raw, y_raw, dimL, want):
    '''
    Function obtain type, x and y coordinate of the desired type:
        type 1 for hard sphere, type 2 for interface 1 and type 3 for interface 2
    
    '''
    
    modType = typ[np.where(typ == want)]
    x = x_raw[np.where(typ == want)] 
    y = y_raw[np.where(typ == want)]
    
    return modType,x, y




def checkAngle(atoms, boxL, xSticky2, ySticky2, xCenter, yCenter, xSticky1, ySticky1):
    for atom in range(atoms):

        A = (Dist(xSticky2[atom],xCenter[atom], boxL)**2 + Dist(ySticky2[atom],yCenter[atom], boxL)**2)**(1/2)
        B = (Dist(xSticky1[atom],xCenter[atom], boxL)**2 + Dist(ySticky1[atom],yCenter[atom], boxL)**2)**(1/2)
        C = (Dist(xSticky1[atom],xSticky2[atom], boxL)**2 + Dist(ySticky1[atom],ySticky2[atom], boxL)**2)**(1/2)

        numerator = A**2 + B**2 - C**2
        denom = 2 * A * B

        angle = np.arccos(numerator/denom)
        print(angle)
    


def RDF_Snap(tag, dimL, edges):
    dist = []
    for atom in range(atoms):
        for atom2 in range(atoms):
            if atom != atom2:

                dx = Dist(tag[atom][0],tag[atom2][0], dimL)
                dy = Dist(tag[atom][1], tag[atom2][1], dimL)
                               
                dist.append((dx**2 + dy**2)**0.5)

    rdf, edges = np.histogram(dist, bins = edges, density = False)

    return rdf
def combine_array(x,y):
    array = []
    
    for i in range(len(x)):
        array.append([x[i], y[i]])
    return np.array(array)

def RDF_SnapCenters(x, y, box1D, dimL, edges):
    dist = []
    for atom in range(atoms):
        for atom2 in range(atoms):
            if atom != atom2:

                dx = Dist(x[atom], x[atom2], dimL)
                dy = Dist(y[atom], y[atom2], dimL)
                                   
                dist.append((dx**2 + dy**2)**0.5)

    rdf, edges = np.histogram(dist, bins = edges, density = False)

    return rdf
def CalcLJ(sigma, r): #sigma is cut off distance
    atomLJ = ((sigma/r)**12) - ((sigma/r)**6)   
        
    return atomLJ


def findtagangle(tag,boxlength):
    center =[0,0]
    point2 = [4,0]
    tag = np.array(tag, dtype = 'float')
    
   
    costag_all = []    
    for i in range(len(tag)):
        
        #Between tag - center ( 0,0)
        xdis = Periodic(center[0]- tag[i][0],boxlength)
        ydis = Periodic(center[1] - tag[i][1],boxlength)
        
        AC = (xdis**2 + ydis**2 )**0.5
        #between 2 point in x axis
        xdis = Periodic(center[0]- point2[0],boxlength)
        ydis = Periodic(center[1] - point2[1],boxlength)
        
        AB = (xdis**2 + ydis**2) **0.5
        #between tag- random chosen point on x axis
        xdis = Periodic(point2[0] - tag[i][0],boxlength)
        ydis = Periodic(point2[1] - tag[i][1],boxlength)
        
        BC = (xdis **2 + ydis**2) **0.5
        
        num = BC**2 -AC**2 -AB**2
        den = -2 *AC * AB
        costag = num/den
        costag_all.append(costag)
    return costag_all
      
#simulation info
'''
Have to input:
    
    Atom numbers 
    File name to read ( in .txt format)
    Total Move (propMove): Must be equal or less than total moves used in the simulation
    Size of window ( Window) : Must be less than "allsteps = propMove/dump"
    



'''

atoms = 169
file = 'LJMonte_0.28274-eps220.00_eps335.00_eps230.00_beta1.000-2.txt'
atoms_withtag = atoms * 3
dump = 100
propMove = 100000
snaps = int(propMove/dump)
boxL = 30

head = 9
interval = 1 #if you want to skip certain snapshots 
box1D = 4 #number of cells in 1D

#lj parameters 
sigma = 0.25
cutoff = 2.5
epsTag22 = 0.5
epsTag33 = 0.5
epsTag23 = 1.5
diameter_1 = 1 #5 nm
diameterTag2 = 0.2
diameterTag3 = 0.2

#rdf conditions
usedSnap = 1
avgSnap = int(usedSnap * snaps) #double check

bins = 85
maxR = 15 #5 * diameter_1

readSnaps = ReadSnaps(file, head, atoms_withtag, snaps, interval)
startSnap = snaps - avgSnap
val = 0
r_center_int2 = 0.5


Angle_AllStep = []

for snap in range(snaps): #calculate eachsnap

    mol, typ, x, y = next(readSnaps)

    
    if snap in range(startSnap, snaps):
        #print(len(Type_DataCenter(typ,x,y,boxL,1)))
                                                                                                                                                  
        
        
        Type2, xSticky2, ySticky2 = Type_Data(typ, x, y, boxL, 2)
        
        
        sticky2_coord = combine_array(xSticky2, ySticky2)
        angle_eachstep = findtagangle(sticky2_coord,boxL)
        Angle_AllStep.append(angle_eachstep)
        
        
window = int(0.5*propMove/dump)
allsteps = int(propMove/dump)
correlationfunction = []
track = 0
for t_0 in range(allsteps):    
    
    angle_time0 = Angle_AllStep[t_0]
    angle_eachwindow = []
    if t_0 + window >= (allsteps):
        print(t_0)
        break
    else:
        for t in range(window): #for each point in window
            angle_timet = np.array(Angle_AllStep[t_0+t]) # time i
            avgangle_eachpoint = 0
            for eachparticle in range(len(angle_time0)):
                avgangle_eachpoint += angle_time0[eachparticle] * angle_timet[eachparticle]
            angle_eachwindow.append(avgangle_eachpoint/len(angle_time0))
    track +=1 
    print(track)
    correlationfunction.append(angle_eachwindow)
    
final_cfunction = np.mean(correlationfunction,axis = 0)
t = np.arange(0,window*dump,dump)
                
plt.figure()
plt.plot(t,final_cfunction)
plt.xlabel('Time(Steps)')
plt.ylabel('<cos$\Theta$(0)*cos$\Theta$(t)>')
plt.yticks([0,0.1,0.2,0.3,0.4,0.5])

plt.show()


       
def cal_dis(x2,x1,y2,y1):
    return ((x2 - x1)**2 + (y2 - y1)**2 )**0.5
