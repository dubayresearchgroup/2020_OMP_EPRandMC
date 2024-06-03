import random
import math
import numpy as np

def Distance(v1,v2, boxlength):
    dist = abs(v1-v2)
    if dist > (0.5 * boxlength):
        dist = boxlength - dist

    return dist

def CalcLJ(sigma, r): #sigma is diameter of interface
    atomLJ = ((sigma/r)**12) - ((sigma/r)**6)   
        
    return atomLJ

def CreateParticle(boxlength, no_of_atom, diameter, diameterTag2, diameterTag3, sigma, eps22, eps33, eps23, cutoff, angle):
    '''
    To create an evenly spaced array of particle positions. The no_of_atom can only be the number that can be squared
    '''
    
    x_possible = np.zeros(int(np.sqrt(no_of_atom)))
    y_possible = np.zeros(int(np.sqrt(no_of_atom)))
    
    centralPosition= np.arange(2).reshape(1,2)
    
    for i in range(len(x_possible)):
        
        if i == 0:
            x_possible[i] = 1.5*diameter
            y_possible[i] = 1.5*diameter
        else:
            x_possible[i] += x_possible[i-1] + boxlength/np.sqrt(no_of_atom)
            y_possible[i] += x_possible[i-1] + boxlength/np.sqrt(no_of_atom)
    
    for x in x_possible:
        for y in y_possible:
            singleAtom = np.zeros(2)
            singleAtom[0] = x
            singleAtom[1] = y
            singleAtom = singleAtom.reshape(1,2)
            centralPosition = np.concatenate((centralPosition,singleAtom), axis = 0)
            #creates an array with dimensions of number atoms and 2.
            # call this as centralPosition[atomNum][0] for x and centralPosition[atomNum][1] for 
    centralPosition = np.delete(centralPosition, 0, 0)

    #Getting ready to hold tags
    tag2pos = np.array(centralPosition)
    tag3pos = np.array(centralPosition) 

    #determine new location of tags
    for i in range(no_of_atom):
        
        #theta = random.choice(degree)
        rotate = random.uniform(-180,180)
        rad = (np.pi * rotate)/180
        tag2pos[i][0] += np.cos(rad)*(0.5*diameter- (0.5 * 2**(1/6)*sigma))
        tag2pos[i][1] += np.sin(rad)*(0.5*diameter - (0.5 * 2**(1/6)*sigma)) #sticky 2
        
        theta2 = rotate - angle
        rad2= (np.pi * theta2)/180
        tag3pos[i][0] += np.cos(rad2)*(0.5*diameter- (0.5 * 2**(1/6)*sigma)) #sticky 1s
        tag3pos[i][1] += np.sin(rad2)*(0.5*diameter- (0.5 * 2**(1/6)*sigma))


    #Calculate initial energy
    starting_lj = 0

    for start in range(no_of_atom):
        for second in range(start+1,no_of_atom):

            xdist = Distance(centralPosition[start][0],centralPosition[second][0], boxlength)
            ydist = Distance(centralPosition[start][1],centralPosition[second][1],boxlength)
            
            rij = (xdist*xdist + ydist*ydist) **0.5
            if rij  < diameter:
                print('distance between hard spheres is less than diameter_HS')
            
            xdist = Distance(tag2pos[start][0],tag2pos[second][0], boxlength)
            ydist = Distance(tag2pos[start][1],tag2pos[second][1],boxlength)
            
            rij2 = (xdist*xdist + ydist*ydist) **0.5
            if rij2  <= cutoff:
                starting_lj += 4 * eps22 * CalcLJ(sigma, rij2)

            xdist = Distance(tag3pos[start][0],tag3pos[second][0], boxlength)
            ydist = Distance(tag3pos[start][1],tag3pos[second][1], boxlength)
            
            rij3 = (xdist*xdist + ydist*ydist) **0.5
            if rij3  <= cutoff:
                starting_lj += 4 * eps33 * CalcLJ(sigma, rij3)
                
    for tag2 in range(no_of_atom):
        for tag3 in range(tag2+1, no_of_atom):
            if tag2 != tag3:

                    
                xdist = Distance(tag2pos[tag2][0],tag3pos[tag3][0], boxlength)
                ydist = Distance(tag2pos[tag2][1],tag3pos[tag3][1],boxlength)

                
                rij23 = (xdist*xdist + ydist*ydist) **0.5

                if rij23 <= cutoff:
                    starting_lj += 4 * eps23 * CalcLJ(sigma, rij23)

                xdist = Distance(tag3pos[tag2][0],tag2pos[tag3][0], boxlength)
                ydist = Distance(tag3pos[tag2][1],tag2pos[tag3][1],boxlength)

                
                rij32 = (xdist*xdist + ydist*ydist) **0.5

                if rij32 <= cutoff:
                    starting_lj += 4 * eps23 * CalcLJ(sigma, rij32)

    return centralPosition, tag2pos, tag3pos, starting_lj
