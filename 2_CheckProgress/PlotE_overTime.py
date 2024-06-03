import numpy as np
import matplotlib.pyplot as plt

VolFrac = [0.28274] #0.38485 14748
eps11 = [0, 10, 0, 0, 10, 10]
eps22 = [0, 0, 10, 0, 10, 10]
eps12 = [0, 0, 0, 10, 0, 10] 
beta = 1.000
partL = 3000000

energyCheck = 10000

prints = int(partL / energyCheck)

totEnergy = np.zeros((len(eps22),prints))
TransAccept = np.zeros((len(eps22),prints))
RotAccept = np.zeros((len(eps22),prints))
propMove = np.linspace(0, partL, prints)

for vol in range(len(VolFrac)):
    for eps in range(len(eps22)-1):
        file = '%.5f/Trial1/LJMonte_%.5f-eps22%.2f_eps33%.2f_eps23%.2f_beta%.3f_energy.txt' %(VolFrac[0], VolFrac[0], eps22[eps+1], eps11[eps+1], eps12[eps+1], beta)

        RotAccept[eps][:] =(100* np.genfromtxt(file, skip_header = 1, usecols = 1, max_rows = prints) / np.genfromtxt(file, skip_header = 1, usecols =2, max_rows = prints))

        TransAccept[eps][:] =(100* np.genfromtxt(file, skip_header = 1, usecols = 3, max_rows = prints) / np.genfromtxt(file, skip_header = 1, usecols = 4, max_rows = prints))

        
        totEnergy[eps][:] = np.genfromtxt(file, skip_header = 1, usecols = 5, max_rows = prints)

x = 1
plt.plot(propMove[::x],totEnergy[0][::x], color = 'k', label = '$\\epsilon_{11}$ = %.2f, $\\epsilon_{22}$ = %.2f, $\\epsilon_{12}$= %.2f' %(eps11[0], eps22[0], eps12[0]))
plt.plot(propMove[::x],totEnergy[1][::x], color = '#00CC99', label = '$\\epsilon_{11}$ = %.2f, $\\epsilon_{22}$ = %.2f, $\\epsilon_{12}$= %.2f' %(eps11[1], eps22[1], eps12[1]))
plt.plot(propMove[::x],totEnergy[2][::x], color = '#FF0000', label = '$\\epsilon_{11}$ = %.2f, $\\epsilon_{22}$ = %.2f, $\\epsilon_{12}$= %.2f' %(eps11[2], eps22[2], eps12[2]))
plt.plot(propMove[::x],totEnergy[3][::x], color = '#9900CC', label = '$\\epsilon_{11}$ = %.2f, $\\epsilon_{22}$ = %.2f, $\\epsilon_{12}$= %.2f' %(eps11[3], eps22[3], eps12[3]))
plt.plot(propMove[::x],totEnergy[4][::x], color = '#FF6600', label = '$\\epsilon_{11}$ = %.2f, $\\epsilon_{22}$ = %.2f, $\\epsilon_{12}$= %.2f' %(eps11[4], eps22[4], eps12[4]))
plt.plot(propMove[::x],totEnergy[4][::x], color = '#009900', label = '$\\epsilon_{11}$ = %.2f, $\\epsilon_{22}$ = %.2f, $\\epsilon_{12}$= %.2f' %(eps11[5], eps22[5], eps12[5]))

plt.xticks([0, 1000000, 2000000, 3000000])
#plt.xticks()
#plt.ylim(bottom = -1.75)
#plt.ylim(top = 100)
plt.ylim(top = 0)
plt.xlim(left = 0)
plt.xlim(right = partL)
plt.legend(title='Well Depth ($k_BT$)')
plt.xlabel('Proposed Move')
plt.ylabel('Potential Energy')
plt.show()

        



