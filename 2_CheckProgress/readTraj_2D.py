#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 09:53:34 2019

@author: katerilab
"""

import numpy as np


def ReadLines (filename):

    with open(filename, 'r') as file:

        for line in file:

            yield line

def ReadSnap (lines, atoms, head):

    mol = np.empty(atoms)
    typ = np.empty(atoms)
    x = np.empty(atoms)
    y = np.empty(atoms)

    for line in range(head):

        next(lines)

    for atom in range(atoms):

        line = next(lines)
        #if len(line) != 6:
        #print(line)

        atomNum, molID, typS, xS, yS, zS = line.split()
        mol[atom] = int(molID)
        typ[atom] = int(typS)
        x[atom] = float(xS)
        y[atom] = float(yS)

    return mol, typ, x, y

def ReadSnaps (file, head, atoms, snaps, interval):

    skip = (interval - 1) * (atoms + head)

    lines = ReadLines(file)
 

    snap = 0
    while snap < snaps:
        
    
       
        yield ReadSnap(lines, atoms, head)
        
    
        for line in range(skip):
            next(lines)
       
            
        snap += 1
        
        
        
