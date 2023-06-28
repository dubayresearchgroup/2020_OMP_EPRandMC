import numpy as np
from numpy import pi
import matplotlib.patches as patches
from matplotlib.patches import Patch
import time
import matplotlib.pyplot as plt
from readTraj_2D import ReadSnaps
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar, AnchoredDrawingArea
import matplotlib.font_manager as fm
fontprops = fm.FontProperties(size=8)


def getTypeData(typ, x_raw, y_raw, dimL, want):
    
    modType = typ[np.where(typ == want)]
    x = x_raw[np.where(typ == want)] 
    y = y_raw[np.where(typ == want)]
    
    return modType, x, y

def Dist(p1, p2, boxL):

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
        


#simulation info

boxL = 30
propMove = 3000000
dump = 1000
snaps = int(propMove/dump)


head = 9
interval = 1 #if you want to skip certain snapshots 
box1D = 4 #number of cells in 1D

beta = 1.0 
eps11 = [0]
eps22 = [0]
eps21 = [0] 
atoms = [324]
trials = [1]

avgSnap = int(1 * snaps)

##These are the distances and angles of each tag from interface 2
#order for each array: 65, 90, 93, 188, 237, 451, 491, 534
tagShift = [1.1861097242924719, -2.0904901272346796, -2.5819587472789434, 0.2129156857672716, -0.3718925981451964, -3.0050880038014025, 2.2567852620217157, 1.2782621710148958]
tagDist = [0.0919157, 0.165364, 0.1274755, 0.325724, 0.317077, 0.2221824, 0.416997, 0.2085426]


for vol in range(len(atoms)):
    for poss in range(len(eps21)):
        for tri in range(len(trials)):

            volFrac = (atoms[vol] * pi * (0.5)**2) / (boxL**2)
            totatoms = atoms[vol] * 3

            name = 'LJMonte_%.5f-eps22%.2f_eps33%.2f_eps23%.2f_beta%.3f' %(volFrac, eps22[poss], eps11[poss], eps21[poss], beta)
            file = 'OldTest/%s-2.txt' %(name)
            trajFile = 'trial%i-%s-tag.lammpstrj' %(trials[tri], name)

            readSnaps = ReadSnaps(file, head, totatoms, snaps, interval)
            startSnap = snaps - avgSnap
            val = 0

            test = 0

            for snap in range(snaps):
                dist = 0

                mol, typ, x, y = next(readSnaps)

                if snap in range(startSnap, snaps):
                    #get position of hard sphere and each interface 

                    centerType, xCenter, yCenter = getTypeData( typ, x, y, boxL, 1)
                    Type2, xSticky2, ySticky2 = getTypeData(typ, x, y, boxL, 2)
                    Type1, xSticky1, ySticky1 = getTypeData(typ, x, y, boxL, 3)

                    #get position of each tag
                    tagX65, tagY65 = gettag_coordinateNew(atoms[0], tagShift[0], tagDist[0],  xCenter, yCenter,xSticky2, ySticky2, boxL)
                    tagX90, tagY90 = gettag_coordinateNew(atoms[0], tagShift[1], tagDist[1],  xCenter, yCenter,xSticky2, ySticky2, boxL)
                    tagX93, tagY93 = gettag_coordinateNew(atoms[0], tagShift[2], tagDist[2],  xCenter, yCenter,xSticky2, ySticky2, boxL)
                    tagX188, tagY188 = gettag_coordinateNew(atoms[0], tagShift[3], tagDist[3], xCenter, yCenter,xSticky2, ySticky2, boxL)
                    tagX237, tagY237 = gettag_coordinateNew(atoms[0],tagShift[4], tagDist[4], xCenter, yCenter,xSticky2, ySticky2, boxL)
                    tagX451, tagY451 = gettag_coordinateNew(atoms[0],tagShift[5], tagDist[5], xCenter, yCenter,xSticky2, ySticky2, boxL)
                    tagX491, tagY491 = gettag_coordinateNew(atoms[0],tagShift[6], tagDist[6], xCenter, yCenter,xSticky2, ySticky2, boxL)
                    tagX534, tagY534 = gettag_coordinateNew(atoms[0], tagShift[7], tagDist[7], xCenter, yCenter,xSticky2, ySticky2, boxL)
                    

        fig, ax = plt.subplots()
        ax.set_xlim((0,30))
        ax.set_ylim((0,30))

        for i in range(atoms[0]):
            #plot the center atom
            circ = patches.Circle((xCenter[i],yCenter[i]), 0.5,  facecolor = 'none', edgecolor = 'k', label = 'Center')
            
            #plot the interfaces 
            s2 = patches.Circle((xSticky2[i],ySticky2[i]), 0.2,  facecolor = 'none', edgecolor = 'b', label = 'Sticky 2')
            s1 = patches.Circle((xSticky1[i],ySticky1[i]), 0.2,  facecolor = 'none', edgecolor = 'g', label = 'Sticky 1')
            
            #plot each tag
            tag188 = patches.Circle((tagX188[i],tagY188[i]), 0.1,  facecolor = '#6600cc', edgecolor = '#6600cc', label = 'Tag 188')
            tag534 = patches.Circle((tagX534[i],tagY534[i]), 0.1,  facecolor = '#ff6600', edgecolor = '#ff6600', label = 'Tag 534')
            tag90 = patches.Circle((tagX90[i],tagY90[i]), 0.1,  facecolor = '#ffcc00', edgecolor = '#ffcc00', label = 'Tag 90')
            tag237 = patches.Circle((tagX237[i],tagY237[i]), 0.1,  facecolor = '#339966', edgecolor = '#339966', label = 'Tag 237')
            tag65 = patches.Circle((tagX65[i],tagY65[i]), 0.1,  facecolor = '#ff00ff', edgecolor = '#ff00ff', label = 'Tag 65')
            tag93 = patches.Circle((tagX93[i],tagY93[i]), 0.1,  facecolor = '#00ccff', edgecolor = '#00ccff', label = 'Tag 93')
            tag451 = patches.Circle((tagX451[i],tagY451[i]), 0.1,  facecolor = '#99FF99', edgecolor = '#99FF99', label = 'Tag 451')
            tag491 = patches.Circle((tagX491[i],tagY491[i]), 0.1,  facecolor = '#CC0000', edgecolor = '#CC0000', label = 'Tag 491')
            
            #add each component to the image 
            ax.add_artist(circ)
            ax.add_artist(s2)
            ax.add_artist(s1)
            ax.add_artist(tag188)
            ax.add_artist(tag534)
            ax.add_artist(tag90)
            ax.add_artist(tag237)
            ax.add_artist(tag65)
            ax.add_artist(tag93)
            ax.add_artist(tag451)
            ax.add_artist(tag491)
            
            
        legend_elements = [Line2D([0], [0], marker='o', color='none', label='Center', markerfacecolor='none',markeredgecolor='k', markersize= 15),
                           Line2D([0], [0], marker='o', color='none', label='Sticky 2', markerfacecolor='none',markeredgecolor='b', markersize= 6),
                           Line2D([0], [0], marker='o', color='none', label='Sticky 1', markerfacecolor='none',markeredgecolor='g', markersize= 6),
                           Line2D([0], [0], marker='o', color='none', label='Tag 188', markerfacecolor='#6600cc',markeredgecolor='#6600cc', markersize= 3),
                           Line2D([0], [0], marker='o', color='none', label='Tag 534', markerfacecolor='#ff6600',markeredgecolor='#ff6600', markersize= 3),
                           Line2D([0], [0], marker='o', color='none', label='Tag 90', markerfacecolor='#ffcc00',markeredgecolor='#ffcc00', markersize= 3),
                           Line2D([0], [0], marker='o', color='none', label='Tag 93', markerfacecolor='#00ccff',markeredgecolor='#00ccff', markersize= 3),
                           Line2D([0], [0], marker='o', color='none', label='Tag 65', markerfacecolor='#ff00ff',markeredgecolor='#ff00ff', markersize= 3),
                           Line2D([0], [0], marker='o', color='none', label='Tag 451', markerfacecolor='#99FF99',markeredgecolor='#99FF99', markersize= 3),
                           Line2D([0], [0], marker='o', color='none', label='Tag 491', markerfacecolor='#CC0000',markeredgecolor='#CC0000', markersize= 3),
                           Line2D([0], [0], marker='o', color='none', label='Tag 237', markerfacecolor='#339966',markeredgecolor='#339966', markersize= 3)]
                           #Line2D([0], [0], color='k', lw=1, label='1.1 $\sigma$'), Line2D([0], [0], color='purple', lw=1, label='1.1 $\sigma$')]

        ax.legend(handles=legend_elements)

        scalebar = AnchoredSizeBar(ax.transData,
                                   1, '1 $\sigma$', 'lower left', 
                                   pad=0.4,
                                   color='black',
                                   frameon=True,
                                   size_vertical=0.4,
                                   fontproperties=fontprops)
        ax.add_artist(scalebar)

        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_snap('true')
        ax.set_aspect('equal')
        ax.set_title('$\\epsilon_{11}$ = %.2f $k_BT$ , $\\epsilon_{22}$ = %.2f $k_BT$, $\\epsilon_{12}$= %.2f $k_BT$' %(eps11[poss], eps22[poss], eps21[poss]), color = 'k')
        fig.savefig('%s_Snap.png' %(name))
        plt.show() 

                             
