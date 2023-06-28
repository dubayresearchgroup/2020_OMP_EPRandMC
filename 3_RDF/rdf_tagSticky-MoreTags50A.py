import numpy as np
from numpy import pi
import time
from readTraj_2D import ReadSnaps

def getTypeData(typ, x_raw, y_raw, dimL, want):
    ##return numpy array of just data of specific type
    modType = typ[np.where(typ == want)]
    x = x_raw[np.where(typ == want)] 
    y = y_raw[np.where(typ == want)]
    
    return modType, x, y

def Dist(p1, p2, boxL):
    #calculate distance between an ij pair 
    dist = abs(p1 - p2)
    if dist > 0.5 * boxL:
        dist = boxL - dist

    return dist


def gettag_coordinateNew(atoms, shift, tagD, xCenter, yCenter,xSticky2, ySticky2, boxlength):
    #Places tag based off of the distances provided Oct 2nd and the angle at which sticky2 is
    #rotated. Adds/subtracts difference for each tag from the angle of Sticky2 
    tagX = np.zeros(atoms)
    tagY = np.zeros(atoms)
    for atom in range(atoms):
        
        sticky2Angle = np.arctan2(ySticky2[atom] - yCenter[atom], xSticky2[atom] - xCenter[atom])

        tagAngle = sticky2Angle + shift
        
        tagX[atom] = (xCenter[atom] + np.cos(tagAngle)*(tagD))
        tagY[atom] = (yCenter[atom]  + np.sin(tagAngle)*(tagD))

        dx = Dist(tagX[atom], xCenter[atom], boxL)
        dy = Dist(tagY[atom], yCenter[atom], boxL)

        dist = (dx**2 + dy**2)**(1/2)      

    return tagX, tagY
        

        
def CalcRDF(atoms, x, y, box1D, dimL, edges):
    dist = []
    for atom in range(atoms):
        for atom2 in range(atoms):
            if atom != atom2:

                dx = Dist(x[atom], x[atom2], dimL)
                dy = Dist(y[atom], y[atom2], dimL)
                                   
                dist.append((dx**2 + dy**2)**0.5)

    rdf, edges = np.histogram(dist, bins = edges, density = False)

    return rdf



#simulation info

boxL = 30
propMove = 3_000_000 ##number of moves 
dump = 1000 ##print frequency, line 342 in sim. code
snaps = int(propMove/dump) ##number of snapshots from simulation 


head = 9
interval = 1 #if you want to skip certain snapshots 
box1D = 4 #number of cells in 1D

beta = 1.0 
eps11 = [10]
eps22 = [10]
eps21 = [10]
atoms = [324]
trials = [1] 


#rdf conditions
avgSnap = int(0.833 * snaps) ###take average over last 83.3% of simulation 

bins = 300
maxR = 15 #5 * diameter_1
edges = np.linspace(0, maxR, bins + 1)
binW = edges[1]-edges[0] #bin width
binsC = np.linspace(0.5*binW, maxR - 0.5*binW, bins) #center of each bin


##These are the distances and angles of each tag from interface 2
#order for each array: 65, 90, 93, 188, 237, 451, 491, 534
tagShift = [1.1861097242924719, -2.0904901272346796, -2.5819587472789434, 0.2129156857672716, -0.3718925981451964, -3.0050880038014025, 2.2567852620217157, 1.2782621710148958]
tagDist = [0.0919157, 0.165364, 0.1274755, 0.325724, 0.317077, 0.2221824, 0.416997, 0.2085426]

for vol in range(len(atoms)):
    for poss in range(len(eps22)):
        for tri in range(len(trials)):
            print(eps22[poss], eps11[poss], eps21[poss])
                        
            ##Arrays to store the RDF in for each tag 
            rdf_centers = np.zeros(len(edges) - 1)
            rdf_tag188 = np.zeros(len(edges) - 1)
            rdf_tag534 = np.zeros(len(edges) - 1)
            rdf_tag90 = np.zeros(len(edges) - 1)
            rdf_tag237 = np.zeros(len(edges) - 1)
            rdf_tag65 = np.zeros(len(edges) - 1)
            rdf_tag93 = np.zeros(len(edges) - 1)
            rdf_tag451 = np.zeros(len(edges) - 1)
            rdf_tag491 = np.zeros(len(edges) - 1)

            normRDF_centers = np.zeros(len(edges)-1)
            normRDF_tag188 = np.zeros(len(edges)-1)
            normRDF_tag534 = np.zeros(len(edges)-1)
            normRDF_tag90 = np.zeros(len(edges)-1)
            normRDF_tag237 = np.zeros(len(edges)-1)
            normRDF_tag65 = np.zeros(len(edges)-1)
            normRDF_tag93 = np.zeros(len(edges)-1)
            normRDF_tag451 = np.zeros(len(edges)-1)
            normRDF_tag491 = np.zeros(len(edges)-1)


            volFrac = (atoms[vol] * pi * (0.5)**2) / (boxL**2)
            totatoms = atoms[vol] * 3

            name = 'LJMonte_%.5f-eps22%.2f_eps33%.2f_eps23%.2f_beta%.3f-2' %(volFrac, eps22[poss], eps11[poss], eps21[poss], beta)
            file = '%.5f/trial%i/%s.txt' %( volFrac, trials[tri], name)
            trajFile = 'trial%i-%s-tag.lammpstrj' %(trials[tri], name)
            name2 = '%.5f-eps22%.2f_eps11%.2f_eps21%.2f_beta%.3f' %(volFrac, eps22[poss], eps11[poss], eps21[poss], beta)

            readSnaps = ReadSnaps(file, head, totatoms, snaps, interval)
            startSnap = snaps - avgSnap
            
            val = 0

            test = 0

            for snap in range(snaps):
                dist = 0

                mol, typ, x, y = next(readSnaps)

                if snap in range(startSnap, snaps):
                    
                    #Get positions of interface 1, 2 and the center of the HS 
                    centerType, xCenter, yCenter = getTypeData( typ, x, y, boxL, 1)
                    Type2, xSticky2, ySticky2 = getTypeData(typ, x, y, boxL, 2)
                    Type1, xSticky1, ySticky1 = getTypeData(typ, x, y, boxL, 3)
                    
                    # Find tag locations for new position for each snap 
                    tagX65, tagY65 = gettag_coordinateNew(atoms[0], tagShift[0], tagDist[0],  xCenter, yCenter,xSticky2, ySticky2, boxL)
                    tagX90, tagY90 = gettag_coordinateNew(atoms[0], tagShift[1], tagDist[1],  xCenter, yCenter,xSticky2, ySticky2, boxL)
                    tagX93, tagY93 = gettag_coordinateNew(atoms[0], tagShift[2], tagDist[2],  xCenter, yCenter,xSticky2, ySticky2, boxL)
                    tagX188, tagY188 = gettag_coordinateNew(atoms[0], tagShift[3], tagDist[3], xCenter, yCenter,xSticky2, ySticky2, boxL)
                    tagX237, tagY237 = gettag_coordinateNew(atoms[0],tagShift[4], tagDist[4], xCenter, yCenter,xSticky2, ySticky2, boxL)
                    tagX451, tagY451 = gettag_coordinateNew(atoms[0],tagShift[5], tagDist[5], xCenter, yCenter,xSticky2, ySticky2, boxL)
                    tagX491, tagY491 = gettag_coordinateNew(atoms[0],tagShift[6], tagDist[6], xCenter, yCenter,xSticky2, ySticky2, boxL)
                    tagX534, tagY534 = gettag_coordinateNew(atoms[0], tagShift[7], tagDist[7], xCenter, yCenter,xSticky2, ySticky2, boxL)
                    
                    


                    #Calculate the RDF for the new tag position
                    rdf_centers += CalcRDF(atoms[vol], xCenter, yCenter, box1D, boxL, edges)
                    rdf_tag65 += CalcRDF(atoms[vol], tagX65, tagY65, box1D, boxL, edges)
                    rdf_tag90 += CalcRDF(atoms[vol], tagX90, tagY90, box1D, boxL, edges)
                    rdf_tag93 += CalcRDF(atoms[vol], tagX93, tagY93, box1D, boxL, edges)
                    rdf_tag188 += CalcRDF(atoms[vol], tagX188, tagY188, box1D, boxL, edges)
                    rdf_tag237 += CalcRDF(atoms[vol], tagX237, tagY237, box1D, boxL, edges)
                    rdf_tag451 += CalcRDF(atoms[vol], tagX451, tagY451, box1D, boxL, edges)
                    rdf_tag491 += CalcRDF(atoms[vol], tagX491, tagY491, box1D, boxL, edges)
                    rdf_tag534 += CalcRDF(atoms[vol], tagX534, tagY534, box1D, boxL, edges)
                    
                    

                    val +=1
                    
            print(val)

            for i in range(bins):

                #Normalize the RDFs         
                normRDF_centers[i] = ((rdf_centers[i]/val) * boxL**2) / (2 * pi * edges[i+1] * binW * atoms[vol] * (atoms[vol] - 1))

                normRDF_tag188[i] = ((rdf_tag188[i]/val) * boxL**2) / (2 * pi * edges[i+1] * binW * atoms[vol] * (atoms[vol] - 1))
                normRDF_tag534[i] = ((rdf_tag534[i]/val) * boxL**2) / (2 * pi * edges[i+1] * binW * atoms[vol] * (atoms[vol] - 1))
                normRDF_tag90[i] = ((rdf_tag90[i]/val) * boxL**2) / (2 * pi * edges[i+1] * binW * atoms[vol] * (atoms[vol] - 1))
                normRDF_tag237[i] = ((rdf_tag237[i]/val) * boxL**2) / (2 * pi * edges[i+1] * binW * atoms[vol] * (atoms[vol] - 1))

                normRDF_tag451[i] = ((rdf_tag451[i]/val) * boxL**2) / (2 * pi * edges[i+1] * binW * atoms[vol] * (atoms[vol] - 1))
                normRDF_tag491[i] = ((rdf_tag491[i]/val) * boxL**2) / (2 * pi * edges[i+1] * binW * atoms[vol] * (atoms[vol] - 1))
                normRDF_tag93[i] = ((rdf_tag93[i]/val) * boxL**2) / (2 * pi * edges[i+1] * binW * atoms[vol] * (atoms[vol] - 1))
                normRDF_tag65[i] = ((rdf_tag65[i]/val) * boxL**2) / (2 * pi * edges[i+1] * binW * atoms[vol] * (atoms[vol] - 1))

            #Write RDFS to textfile 
            ##There are two different writes here -- BinsEdge and BinsC. The BinsEdge plots will be more like a traditional histogram,
            ## which will allow you to see jumps more clearly. The BinsC is polished curves more traditionally 
            ## shown for RDFS 
            #TODO : make this a function and loop through
            with open('TextFiles/BinsEdge/Trial%i_RDF_Tag188_%s.txt' %(trials[tri], name2), 'w') as filehandle:
                for value in range(len(normRDF_tag188)):
                
                    filehandle.write('\n')
                    filehandle.write('%.4f ' %(binsC[value]-(0.5*binW)))
                    filehandle.write('%.6f ' %(normRDF_tag188[value]))
                    filehandle.write('\n')
                    filehandle.write('%.4f ' %(binsC[value]+(0.5*binW)))
                    filehandle.write('%.6f ' %(normRDF_tag188[value]))

            
            with open('TextFiles/BinsEdge/Trial%i_RDF_Tag534_%s.txt' %(trials[tri],name2), 'w') as filehandle:
                for value in range(len(normRDF_tag534)):
                
                    filehandle.write('\n')
                    filehandle.write('%.4f ' %(binsC[value]-(0.5*binW)))
                    filehandle.write('%.6f ' %(normRDF_tag534[value]))
                    filehandle.write('\n')
                    filehandle.write('%.4f ' %(binsC[value]+(0.5*binW)))
                    filehandle.write('%.6f ' %(normRDF_tag534[value]))

            
            with open('TextFiles/BinsEdge/Trial%i_RDF_Tag90_%s.txt' %(trials[tri],name2), 'w') as filehandle:
                for value in range(len(normRDF_tag90)):
                
                    filehandle.write('\n')
                    filehandle.write('%.4f ' %(binsC[value]-(0.5*binW)))
                    filehandle.write('%.6f ' %(normRDF_tag90[value]))
                    filehandle.write('\n')
                    filehandle.write('%.4f ' %(binsC[value]+(0.5*binW)))
                    filehandle.write('%.6f ' %(normRDF_tag90[value]))
                    

            
            with open('TextFiles/BinsEdge/Trial%i_RDF_Tag237_%s.txt' %(trials[tri],name2), 'w') as filehandle:
                for value in range(len(normRDF_tag237)):
                
                    filehandle.write('\n')
                    filehandle.write('%.4f ' %(binsC[value]-(0.5*binW)))
                    filehandle.write('%.6f ' %(normRDF_tag237[value]))
                    filehandle.write('\n')
                    filehandle.write('%.4f ' %(binsC[value]+(0.5*binW)))
                    filehandle.write('%.6f ' %(normRDF_tag237[value]))

                    


            with open('TextFiles/BinsEdge/Trial%i_RDF_Tag65_%s.txt' %(trials[tri], name2), 'w') as filehandle:
                for value in range(len(normRDF_tag65)):
                
                    filehandle.write('\n')
                    filehandle.write('%.4f ' %(binsC[value]-(0.5*binW)))
                    filehandle.write('%.6f ' %(normRDF_tag65[value]))
                    filehandle.write('\n')
                    filehandle.write('%.4f ' %(binsC[value]+(0.5*binW)))
                    filehandle.write('%.6f ' %(normRDF_tag65[value]))

            
            with open('TextFiles/BinsEdge/Trial%i_RDF_Tag93_%s.txt' %(trials[tri],name2), 'w') as filehandle:
                for value in range(len(normRDF_tag93)):
                
                    filehandle.write('\n')
                    filehandle.write('%.4f ' %(binsC[value]-(0.5*binW)))
                    filehandle.write('%.6f ' %(normRDF_tag93[value]))
                    filehandle.write('\n')
                    filehandle.write('%.4f ' %(binsC[value]+(0.5*binW)))
                    filehandle.write('%.6f ' %(normRDF_tag93[value]))

            
            with open('TextFiles/BinsEdge/Trial%i_RDF_Tag451_%s.txt' %(trials[tri],name2), 'w') as filehandle:
                for value in range(len(normRDF_tag451)):
                
                    filehandle.write('\n')
                    filehandle.write('%.4f ' %(binsC[value]-(0.5*binW)))
                    filehandle.write('%.6f ' %(normRDF_tag451[value]))
                    filehandle.write('\n')
                    filehandle.write('%.4f ' %(binsC[value]+(0.5*binW)))
                    filehandle.write('%.6f ' %(normRDF_tag451[value]))
                    

            
            with open('TextFiles/BinsEdge/Trial%i_RDF_Tag491_%s.txt' %(trials[tri],name2), 'w') as filehandle:
                for value in range(len(normRDF_tag491)):
                
                    filehandle.write('\n')
                    filehandle.write('%.4f ' %(binsC[value]-(0.5*binW)))
                    filehandle.write('%.6f ' %(normRDF_tag491[value]))
                    filehandle.write('\n')
                    filehandle.write('%.4f ' %(binsC[value]+(0.5*binW)))
                    filehandle.write('%.6f ' %(normRDF_tag491[value]))



            with open('TextFiles/BinsEdge/Trial%i_RDF_Centers_%s.txt' %(trials[tri],name2), 'w') as filehandle:
                for value in range(len(normRDF_centers)):
                
                    filehandle.write('\n')
                    filehandle.write('%.4f ' %(binsC[value]-(0.5*binW)))
                    filehandle.write('%.6f ' %(normRDF_centers[value]))
                    filehandle.write('\n')
                    filehandle.write('%.4f ' %(binsC[value]+(0.5*binW)))
                    filehandle.write('%.6f ' %(normRDF_centers[value]))

            with open('TextFiles/BinsC/Trial%i_RDF_Tag188_%s.txt' %(trials[tri], name2), 'w') as filehandle:
                for value in range(len(normRDF_tag188)):
                
                    filehandle.write('\n')
                    filehandle.write('%.4f ' %(binsC[value]))
                    filehandle.write('%.6f ' %(normRDF_tag188[value]))

            
            with open('TextFiles/BinsC/Trial%i_RDF_Tag534_%s.txt' %(trials[tri],name2), 'w') as filehandle:
                for value in range(len(normRDF_tag534)):
                
                    filehandle.write('\n')
                    filehandle.write('%.4f ' %(binsC[value]))
                    filehandle.write('%.6f ' %(normRDF_tag534[value]))

            
            with open('TextFiles/BinsC/Trial%i_RDF_Tag90_%s.txt' %(trials[tri],name2), 'w') as filehandle:
                for value in range(len(normRDF_tag90)):
                
                    filehandle.write('\n')
                    filehandle.write('%.4f ' %(binsC[value]))
                    filehandle.write('%.6f ' %(normRDF_tag90[value]))
            
            with open('TextFiles/BinsC/Trial%i_RDF_Tag237_%s.txt' %(trials[tri],name2), 'w') as filehandle:
                for value in range(len(normRDF_tag237)):
                
                    filehandle.write('\n')
                    filehandle.write('%.4f ' %(binsC[value]))
                    filehandle.write('%.6f ' %(normRDF_tag237[value]))


            with open('TextFiles/BinsC/Trial%i_RDF_Tag65_%s.txt' %(trials[tri], name2), 'w') as filehandle:
                for value in range(len(normRDF_tag65)):
                
                    filehandle.write('\n')
                    filehandle.write('%.4f ' %(binsC[value]))
                    filehandle.write('%.6f ' %(normRDF_tag65[value]))

            
            with open('TextFiles/BinsC/Trial%i_RDF_Tag93_%s.txt' %(trials[tri],name2), 'w') as filehandle:
                for value in range(len(normRDF_tag93)):
                
                    filehandle.write('\n')
                    filehandle.write('%.4f ' %(binsC[value]))
                    filehandle.write('%.6f ' %(normRDF_tag93[value]))

            
            with open('TextFiles/BinsC/Trial%i_RDF_Tag451_%s.txt' %(trials[tri],name2), 'w') as filehandle:
                for value in range(len(normRDF_tag451)):
                
                    filehandle.write('\n')
                    filehandle.write('%.4f ' %(binsC[value]))
                    filehandle.write('%.6f ' %(normRDF_tag451[value]))
            
            with open('TextFiles/BinsC/Trial%i_RDF_Tag491_%s.txt' %(trials[tri],name2), 'w') as filehandle:
                for value in range(len(normRDF_tag491)):
                
                    filehandle.write('\n')
                    filehandle.write('%.4f ' %(binsC[value]))
                    filehandle.write('%.6f ' %(normRDF_tag491[value]))

            
            with open('TextFiles/BinsC/Trial%i_RDF_Centers_%s.txt' %(trials[tri],name2), 'w') as filehandle:
                for value in range(len(normRDF_centers)):
                
                    filehandle.write('\n')
                    filehandle.write('%.4f ' %(binsC[value]))
                    filehandle.write('%.6f ' %(normRDF_centers[value]))s
          

