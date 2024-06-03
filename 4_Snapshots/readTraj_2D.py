import numpy as np
##File to iterate through a LAMMPS trajectory 
def ReadLines (filename):

    with open(filename, 'r') as file:

        for line in file:

            yield line

def ReadSnap (lines, atoms, head):

    typ = np.empty(atoms)
    x = np.empty(atoms)
    y = np.empty(atoms)

    for line in range(head):

        next(lines)

    for atom in range(atoms):

        line = next(lines)

        IDstr, typS, xS, yS, zS = line.split()
        typ[atom] = int(typS)
        x[atom] = float(xS)
        y[atom] = float(yS)

    return typ, x, y

def ReadSnaps (file, head, atoms, snaps, interval):

    skip = (interval - 1) * (atoms + head)
    lines = ReadLines(file)

    snap = 0
    while snap < snaps:
        
        yield ReadSnap(lines, atoms, head)
        for line in range(skip):
            next(lines)
            
        snap += 1