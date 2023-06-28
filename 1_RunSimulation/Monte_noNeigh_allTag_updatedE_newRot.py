import random
import math
import numpy as np
from CreateSysAllTag_newRot import CreateParticle
import time
from WriteData import LAMMPS_Initial, analysis_Initial, LAMMPS_step, analysis_step

#! Note : here we utilize tag and interface interchangeably. In this MC code, both 
#! refer to the beta carbons detailed in the publication. Additionally, the interfaces are 
#! labeled numbers 2 and 3. 

## Calculate the LJ MC for a simulation with either attractive tag 2 and 3.  Both contribute to the E, but do not interact.
def Distance(v1,v2, boxL):
    dist = abs(v1-v2)
    if dist > (0.5 * boxL):
        dist = boxL - dist

    return dist

def CalcLJ(sigma, r): #sigma is cut off distance
    atomLJ = ((sigma/r)**12) - ((sigma/r)**6)   
        
    return atomLJ

#%%
def Calc_energyAll(no_of_atom, tag2pos, tag3pos, eps22, eps33, eps23, sigma, boxL):
    ##calaculates the total energy of system for acceptance criteria 
    initial_LJ_all = 0
    for start in range(no_of_atom):
        for second in range(start+1,no_of_atom):
            xdist = Distance(tag2pos[start][0],tag2pos[second][0], boxL)
            ydist = Distance(tag2pos[start][1],tag2pos[second][1],boxL)

            
            rij2 = (xdist*xdist + ydist*ydist) **0.5

            if rij2 < cutoff:
                initial_LJ_all += 4 * eps22 * CalcLJ(sigma, rij2)
                

            xdist = Distance(tag3pos[start][0],tag3pos[second][0], boxL)
            ydist = Distance(tag3pos[start][1],tag3pos[second][1],boxL)
            
            rij3 = (xdist*xdist + ydist*ydist) **0.5

            if rij3 < cutoff:
                initial_LJ_all += 4 * eps33 * CalcLJ(sigma, rij3)
                
            
            xdist = Distance(tag2pos[start][0],tag3pos[second][0], boxL)
            ydist = Distance(tag2pos[start][1],tag3pos[second][1],boxL)

        
            rij23 = (xdist*xdist + ydist*ydist) **0.5

            if rij23 < cutoff:
                initial_LJ_all += 4 * eps23 * CalcLJ(sigma, rij23)
                

            xdist = Distance(tag3pos[start][0],tag2pos[second][0], boxL)
            ydist = Distance(tag3pos[start][1],tag2pos[second][1],boxL)

        
            rij32 = (xdist*xdist + ydist*ydist) **0.5

            if rij32 < cutoff:
                initial_LJ_all += 4 * eps23 * CalcLJ(sigma, rij32)
                
                        
    return initial_LJ_all

#%%
def BackInBox(potentialNewPos, potentialTag2, potentialTag3, movingAtom, boxL):
    ##Wrap particle back into box if it crosses the periodic boundary condition 
    for i in range(2):
        #check x-component 
        if potentialNewPos[movingAtom][i] > boxL:
            potentialNewPos[movingAtom][i] -= boxL
            potentialTag2[movingAtom][i] -= boxL
            potentialTag3[movingAtom][i] -= boxL
        #check y-component 
        if potentialNewPos[movingAtom][i] < 0:
            potentialNewPos[movingAtom][i] += boxL
            potentialTag2[movingAtom][i] += boxL
            potentialTag3[movingAtom][i] += boxL

#%%
def CalcE_moveAtom(no_of_atom, movingAtom, tag2, tag3, eps22, eps33, eps23, sigma, boxL,cutoff):
    ##Calculate the energy from a single move 
    E = 0
    for neigh in range(no_of_atom):
        #TODO : impliment neighborlist for speed up
        if neigh != movingAtom:
            
            #calculate distance for interface #TODO #? 
            dx22 = Distance(tag2[movingAtom][0],tag2[neigh][0], boxL)
            dy22 = Distance(tag2[movingAtom][1],tag2[neigh][1], boxL)
            dist22 = (dx22**2 + dy22**2)**(1/2)

            if dist22 < cutoff:
                E += 4 * eps22 * CalcLJ(sigma, dist22)

            #calculate distance for interface #TODO #? 
            dx33 = Distance(tag3[movingAtom][0],tag3[neigh][0], boxL)
            dy33 = Distance(tag3[movingAtom][1],tag3[neigh][1], boxL)
            dist33 = (dx33**2 + dy33**2)**(1/2)

            if dist33 < cutoff:
                E += 4 * eps33 * CalcLJ(sigma, dist33)

            #calculate distance for interface 1 - 2 
            dx23 = Distance(tag2[movingAtom][0],tag3[neigh][0], boxL)
            dy23 = Distance(tag2[movingAtom][1],tag3[neigh][1], boxL)
            dist23 = (dx23**2 + dy23**2)**(1/2)

            if dist23 < cutoff:
                E += 4 * eps23 * CalcLJ(sigma, dist23)

                ##This section does not double count -- without it
                ## only half of the 1-2 energetics are calculated 
            dx32 = Distance(tag3[movingAtom][0],tag2[neigh][0], boxL)
            dy32 = Distance(tag3[movingAtom][1],tag2[neigh][1], boxL)
            dist32 = (dx32**2 + dy32**2)**(1/2)

            if dist32 < cutoff:
                E += 4 * eps23 * CalcLJ(sigma, dist32)

    return E

#%%           
     
    
def LJMonteMove(centralPosition, tag2pos, tag3pos, boxL, num_of_atom, dump, move,angleMove, diameter_HS, diameter_int1, eps22, eps33, eps23, beta, start, angle,name): 
 
    TranslateAccept = 0
    RotationAccept = 0
    typeR = 0
    typeT = 0
    TotalEnergy = start
    
    optionalStep = ['rotation', 'translate']
    #TODO : break rotation and translate moves into individual functions 

    for step in range(0,steps):
    

        stepType = random.choice(optionalStep)
        movingAtom = random.randint(0,num_of_atom-1)

        initial_LJ_all = start
        deltaE_subsec = 0
                     
        if stepType == 'translate':
            typeT += 1
            deltaEAll = 0
            
            potentialNewPos = np.array(centralPosition)
            potentialTag2 = np.array(tag2pos)
            potentialTag3 = np.array(tag3pos) 

            #determine random translation in x and y directions 
            xshift = random.uniform(-move,move)
            yshift = random.uniform(-move,move)

            #Shift HS and interfaces to new position 
            potentialNewPos[movingAtom][0] += xshift
            potentialNewPos[movingAtom][1] += yshift

            potentialTag2[movingAtom][0] +=  xshift
            potentialTag2[movingAtom][1] +=  yshift
            
            potentialTag3[movingAtom][0] +=  xshift
            potentialTag3[movingAtom][1] +=  yshift

            #check that move did not push particle outside of the box
            #if so, wrap back into box 
            BackInBox(potentialNewPos, potentialTag2, potentialTag3, movingAtom, boxL)
            
            ##Calculate distance 
            for start in range(num_of_atom):
                for second in range(start+1,num_of_atom):
                    xdist = Distance(potentialNewPos[start][0],potentialNewPos[second][0], boxL)
                    ydist = Distance(potentialNewPos[start][1],potentialNewPos[second][1],boxL)

                    rij = (xdist*xdist + ydist*ydist) **0.5

                    if rij < diameter_HS:
                        break
                if rij < diameter_HS:
                    break


            if rij > diameter_HS:
                #if particles are not overlapping
                #calculate the change in energy 

                E_Before = CalcE_moveAtom(num_of_atom, movingAtom, tag2pos, tag3pos, eps22, eps33, eps23, sigma, boxL,cutoff)

                E_After = CalcE_moveAtom(num_of_atom, movingAtom, potentialTag2, potentialTag3, eps22, eps33, eps23, sigma, boxL,cutoff)
                                              

                deltaE_subsec = E_After - E_Before 
                
                if deltaE_subsec <= 0: 
                    #if energy goes down, favorable move is accepted 
                    tag2pos = np.array(potentialTag2)
                    tag3pos = np.array(potentialTag3)
                    
                    TotalEnergy += deltaE_subsec

                    centralPosition = np.array(potentialNewPos)

                    TranslateAccept  +=1

                elif deltaE_subsec > 0:
                    #if energy is  positive, check metropolis 
                    w = np.exp(-beta * deltaE_subsec)
                    randomNum = random.random()

                    if randomNum <= w:
                        tag2pos = np.array(potentialTag2)
                        tag3pos = np.array(potentialTag3)
                        TotalEnergy += deltaE_subsec

                        centralPosition = np.array(potentialNewPos)

                        TranslateAccept  +=1

    
                    
        if stepType == "rotation":
            typeR += 1
            potentialTag2 = np.array(tag2pos)
            potentialTag3 = np.array(tag3pos) 

            ## Position of tags is only updated if r = 0 ##
            theta = random.uniform(-angleMove,angleMove)
            

            ### calculate old theta
            calcAngle = np.arctan2(tag2pos[movingAtom][1] - centralPosition[movingAtom][1], tag2pos[movingAtom][0] - centralPosition[movingAtom][0])

            calcAngle = (calcAngle * 180)/np.pi

            newTotTheta = (calcAngle) + (theta)
            rad = (np.pi * newTotTheta)/180

            #apply rotation to tag 2
            potentialTag2[movingAtom][0] = (centralPosition[movingAtom][0] + np.cos(rad)*(0.5*diameter_HS - (0.5 * 2**(1/6)*sigma)))
            potentialTag2[movingAtom][1] = (centralPosition[movingAtom][1]  +np.sin(rad)*(0.5*diameter_HS - (0.5 * 2**(1/6)*sigma)))

            #calculate rotation on tag 3 and apply it 
            potentialTag3 = np.array(tag3pos)
            theta2 = newTotTheta - angle
            rad2 = (np.pi * theta2)/180

            potentialTag3[movingAtom][0] = (centralPosition[movingAtom][0] + np.cos(rad2)*(0.5*diameter_HS - (0.5 * 2**(1/6)*sigma)))
            potentialTag3[movingAtom][1] = (centralPosition[movingAtom][1] +np.sin(rad2)*(0.5*diameter_HS - (0.5 * 2**(1/6)*sigma)))

            
            ## calculate energ y
           
            E_Before = CalcE_moveAtom(num_of_atom, movingAtom, tag2pos, tag3pos, eps22, eps33, eps23, sigma, boxL,cutoff)

            E_After = CalcE_moveAtom(num_of_atom, movingAtom, potentialTag2, potentialTag3, eps22, eps33, eps23, sigma, boxL,cutoff)
                                              

            deltaE_subsec = E_After - E_Before
        

            if  deltaE_subsec <= 0:
                #if energy goes down, favorable move is accepted 

                tag2pos = np.array(potentialTag2)
                tag3pos = np.array(potentialTag3)
                
                TotalEnergy += deltaE_subsec

                RotationAccept += 1
                
            elif deltaE_subsec > 0:
                #if energy is  positive, check metropolis 

                w = np.exp(-beta * deltaE_subsec)
                randomNum = random.random()

                if randomNum <= w:
                    tag2pos = np.array(potentialTag2)
                    tag3pos = np.array(potentialTag3)
                    TotalEnergy += deltaE_subsec

                    RotationAccept += 1

                
        if (step+1) % (check) == 0 : #print acceptance criteria 
            ReCalcE = Calc_energyAll(num_of_atom, tag2pos, tag3pos, eps22, eps33, eps23, sigma, boxL) #make sure this gets written to something.
            TotalEnergy = ReCalcE

            with open('%s_energy.txt' %(name), 'a') as filehandle:
                filehandle.write('\n')
                filehandle.write('%i %i %i %i %i %.5f' %(step+1, RotationAccept, typeR, TranslateAccept, typeT, TotalEnergy))
           

        if (step+1) % (dump) == 0 : #dump trajectory 
            LAMMPS_step(filename1,num_of_atom,boxL,centralPosition,tag2pos, tag3pos,step)
            analysis_step(filename2,num_of_atom,boxL,centralPosition,tag2pos, tag3pos,step)   
        print(step, stepType)
        
            
    return RotationAccept, TranslateAccept, typeR, typeT


start = time.time()
boxL = 30 #edge of box in units of sigma 
numAtom = [324]
diameter_HS = 1  #diameter of the hardsphere particle
diameter_int1 = 0.20  #diameter of interface 1
diameter_int2 = 0.20 #diameter of interface 2 

## Attraction strength of interface interactions 
## These are epsilon values of the LJ potential 
eps11 = [0, 5.0]
eps22 = [0, 0.0]
eps12 =[0, 0.0]

beta = 1

##This is written in a loop format, so that when one MC simulation finishes the next parameter set would be able to automatically run
for value in range(len(numAtom)):
    for eps in range(len(eps11)):


        #Setting LJ parameter values for current simulation loop 
        a = int(numAtom[value])
        volFrac = (a * math.pi * (0.5 * diameter_HS)**2) / (boxL)**2
        sigma = 0.25
        cutoff = 2.5
        epsTag22 = eps22[eps]  
        epsTag33 = eps11[eps] 
        epsTag23 = eps12[eps]  

        angle = 167 ##Angle between the two interfaces 

        ##Simulation parameters 
        steps = 3000000 ##number of moves 
        dump = 1000 ##print trajectory every this many steps 
        check = 10000 ##print acceptance information every this many steps 
        move = diameter_HS * 0.5 ## maximum length in sigma a particle can move in a single step (moves between [-move,+move])
        angleMove = 45 # maximum angle a particle can rotate in a single move (rotates between [-angleMove,+angleMove])
        

        filename1 = 'LJMonte_%.5f-eps22%.2f_eps33%.2f_eps23%.2f_beta%.3f-2.lammpstrj' %(volFrac, epsTag22, epsTag33, epsTag23, beta)
        filename2 = 'LJMonte_%.5f-eps22%.2f_eps33%.2f_eps23%.2f_beta%.3f-2.txt' %(volFrac, epsTag22, epsTag33, epsTag23, beta)
        name = 'LJMonte_%.5f-eps22%.2f_eps33%.2f_eps23%.2f_beta%.3f' %(volFrac, epsTag22, epsTag33, epsTag23, beta)
        
        ## Get initial particle positions 
        InitialCentralPosition, InitialTag2pos, InitialTag3pos, init_energy = CreateParticle(boxL, a, diameter_HS, diameter_int1, diameter_int2, sigma, epsTag22, epsTag33, epsTag23, cutoff, angle)

        LAMMPS_Initial(filename1, a,boxL,InitialCentralPosition, InitialTag2pos, InitialTag3pos) #this write has no molecule number 
        analysis_Initial(filename2,a,boxL,InitialCentralPosition, InitialTag2pos, InitialTag3pos) #this write has a molecule number 


        ##Run simulation 
        RotationAccept, TranslateAccept, typeR, typeT = LJMonteMove(InitialCentralPosition, InitialTag2pos, InitialTag3pos, boxL, int(a), dump, move,angleMove, diameter_HS, diameter_int1, epsTag22, epsTag33, epsTag23, beta, init_energy, angle, name)
        
        ##Print out acceptance criteria to terminal if on-the-fly tracking is wanted
        print('Rotation Accept = %i, Translation Accept = %i, Total Accept = %i, Percent Accept = %.2f' %(RotationAccept, TranslateAccept, (RotationAccept + TranslateAccept), ((RotationAccept + TranslateAccept)/steps)))
        print('time = %.4f seconds' %(time.time() - start))
        print('Rotation attempt = %i, Translate attempt = %i' %(typeR, typeT))
        print('%.5f %.2f %.2f %.2f %.2f %.2f  %.2f' %(volFrac, epsTag22, epsTag33, epsTag23, (100* (RotationAccept/typeR)), (100 * (TranslateAccept/typeT)), ((RotationAccept + TranslateAccept)/steps)))
        
        ##Save acceptance criteria in a text file 
        with open('AcceptCriteria.txt', 'a') as filehandle:
            filehandle.write('\n')
            filehandle.write('%.5f %.2f %.2f %.2f %.2f %.2f  %.2f' %(volFrac, epsTag22, epsTag33, epsTag23, (100* (RotationAccept/typeR)), (100 * (TranslateAccept/typeT)), (100*(RotationAccept + TranslateAccept)/steps)))



#%%
