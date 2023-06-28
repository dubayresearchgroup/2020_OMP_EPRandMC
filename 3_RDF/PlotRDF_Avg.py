import numpy as np
from numpy import pi
import time
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
fontprops = fm.FontProperties(size=8)

boxL = 30
beta = 1.0 
eps11 = [0, 5, 0, 0, 5, 5]
eps22 = [0, 0, 5, 0, 5, 5]
eps21 = [0, 0, 0,5 , 0, 5]
trials = [1]

#colors for plotting 
allsameInt = ['#00ffcc', '#00b2db', '#0080e6', '#004cf0','#0000ff' ]
likeIntzero = ['#ffccff', '#f08fb2', '#e66680','#db3d4c','#cc0000'] 
unlikeIntzero = ['#ccff33','#99cc26', '#66991a','#33660d','#003300']

atoms = [324]
volFrac = (atoms[0] * pi * (0.5)**2) / (boxL**2)

#tags we are plotting and color 
tags = [188, 534, 90, 237]
tagColor = ['#6600cc','#ff6600','#00CC00','#0066cc' ]


bins = 300
maxR = 15 #5 * diameter_1
edges = np.linspace(0, maxR, bins + 1)
binW = edges[1]-edges[0] #bin width
binsC = np.linspace(0.5*binW, maxR - 0.5*binW, bins) #center of each bin

RDF_int = np.zeros((len(eps11),125))

#get all RDF textfiles and average them together 
for tag in range(len(tags)):
    for value in range(len(eps11)):
        
        trialRDF = np.zeros((1 ,125))
        for tri in range(len(trials)):
          

            name = '%.5f-eps22%.2f_eps11%.2f_eps21%.2f_beta%.3f' %(volFrac, eps22[value], eps11[value], eps21[value], beta)
        
            file = 'Test/Avg_RDF_Tag%i_%.5f_-eps22-%.2f_eps11-%.2f_eps21-%.2f_beta%.3f.txt' %( tags[tag], volFrac, eps22[value], eps11[value], eps21[value], beta)
            bins = np.genfromtxt(file, skip_header = 0, usecols = 0)
        


        RDF_int[value][:] = np.genfromtxt(file, skip_header = 0, usecols = 1) #np.mean(trialRDF, axis=0)
        avg = np.mean(trialRDF, axis=0)

    #plot averaged RDF 
    plt.plot(bins, RDF_int[0], color = 'k', label = '$\\epsilon_{11}$ = %.2f  , $\\epsilon_{22}$ = %.2f , $\\epsilon_{12}$= %.2f' %(eps11[0], eps22[0], eps21[0]))      
    plt.plot(bins, RDF_int[1], color = '#00CC99', label = '$\\epsilon_{11}$ = %.2f , $\\epsilon_{22}$ = %.2f , $\\epsilon_{12}$= %.2f ' %(eps11[1], eps22[1], eps21[1]))
    plt.plot(bins, RDF_int[2], color = '#FF0000',  label = '$\\epsilon_{11}$ = %.2f  , $\\epsilon_{22}$ = %.2f , $\\epsilon_{12}$= %.2f ' %(eps11[2], eps22[2], eps21[2]))
    plt.plot(bins, RDF_int[3], color = '#9900CC',  label = '$\\epsilon_{11}$ = %.2f , $\\epsilon_{22}$ = %.2f , $\\epsilon_{12}$= %.2f ' %(eps11[3], eps22[3], eps21[3]))
    #plt.plot(bins, RDF_int[4], color = '#FF6600', label =  '$\\epsilon_{11}$ = %.2f , $\\epsilon_{22}$ = %.2f , $\\epsilon_{12}$= %.2f ' %(eps11[4], eps22[4], eps21[4]))
    #plt.plot(bins, RDF_int[5], color = '#009900', label = '$\\epsilon_{11}$ = %.2f , $\\epsilon_{22}$ = %.2f , $\\epsilon_{12}$= %.2f ' %(eps11[5], eps22[5], eps21[5]))
    #plt.scatter(binsC[peakTagLoc], avgTag[peakTagLoc], color = 'purple')
    #plt.annotate('%.4f $r/\sigma$ = $\AA$' %(binsC[peakTagLoc], (binsC[peakTagLoc] * 40)), (binsC[peakTagLoc], avgTag[peakTagLoc]))
    plt.xlabel('$r/\sigma$')
    plt.ylabel('g($r\,$)')
    plt.xlim(left = 0)
    plt.xlim(right = 3)
    plt.xticks([0, 0.5, 1, 1.5, 2, 2.5, 3])
    #plt.title('$\\epsilon_{11}$ = %.2f $k_BT$ , $\\epsilon_{22}$ = %.2f $k_BT$, $\\epsilon_{12}$= %.2f $k_BT$' %(eps1122[poss], eps1122[poss], eps21))
    plt.ylim(bottom = 0)

    plt.legend(title='Well Depth ($k_BT$)')
    plt.ylim(top = 4.0)
    title_obj = plt.title('Radial Distribution Function for Tag %i' %(tags[tag]))
    plt.getp(title_obj)                    #print out the properties of title
    plt.getp(title_obj, 'text')            #print out the 'text' property for title
    plt.setp(title_obj, color= '%s' %(tagColor[tag]))
 #  plt.savefig('BinsEdge_tag%i-%.5f-weaker-oneInt_RDF.png' %( tags[tag], volFrac))
    plt.show()



        

        
                   

                
        
