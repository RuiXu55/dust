import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from pylab import *
from matplotlib import rc, font_manager
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable


au = '50'
s1 = '1'
out_ind   = [s1 +'m'+au+'auv.tab']
lab =['-','--',':']
############################################################################


# find the name of the plotting
for l in range(0,1):
  name = 'init_neutral_mass'
  for s in range(0,1):
    if s==0:
      f = open(out_ind[l],'r')
    else:
      f = open(out_ind1[l],'r')
    mark = True
    while(mark):
     line = f.readline()
     data = line.split()
     if data[0]==name:
       mark = False
       numx = int(data[1])
       print ('numx',numx)
       px = np.linspace(0,0,numx)
       neu = np.linspace(0,0,numx)
       for i in range(0,numx):
         line = f.readline()
         data = line.split()
         px[i] = float(data[0])
         neu[i] = float(data[1])
       break
    f.close()

    f = open(out_ind[l],'r')
    numx = numx -1
    name = 'adsorption'
    mark = True
    while(mark):
     line = f.readline()
     data = line.split()
     if data[0]==name:
       mark = False
       col = []
       numy = int(data[2])
       print ('numy',numy)

       adr = np.ones((numx,numy))
       desr = np.ones((numx,numy))
       adf = np.ones((numx,numy))
       desf = np.ones((numx,numy))
       for i in range(0,numx-2):
         fi = f.readline()
         data = fi.split()
         for k in range(0,numy):
             print ('i=',i,' k=',k)
             adr[i,k] = abs(float(data[k]))
         for k in range(0,numy):
             desr[i,k] = abs(float(data[k+numy]))
         for k in range(0,numy):
             adf[i,k] = abs(float(data[k+2*numy]))
         for k in range(0,numy):
             desf[i,k] = abs(float(data[k+3*numy]))

######################## Plot the results #############################
Z = adf.transpose()
print ('z',Z)
extents = (0.0,5,1.0,1e4)
cmap = cm.get_cmap('Greys_r',20)    # 11 discrete colors
cmap = cm.get_cmap('jet')    # 11 discrete colors
#myplt = plt.imshow(Z,extent= extents,cmap=cmap,interpolation="bicubic",
#      aspect='auto',origin='lower',norm=LogNorm(vmin=0.0001, vmax=0.1))
myplt = plt.imshow(Z,extent= extents,cmap=cmap,interpolation="bicubic",
      aspect='auto',origin='lower')

plt.yscale('log')
sizeOfFont = 12
fontProperties = {'family':'sans-serif','sans-serif':['Helvetica'],
        'weight' : 'normal', 'size' : sizeOfFont}
ticks_font = font_manager.FontProperties(family='Helvetica', style='normal',
        size=sizeOfFont, weight='normal', stretch='normal')
rc('text', usetex=True)
rc('font',**fontProperties)
font = {'size'  : 24}
rc('font', **font)
plt.xlabel(r'$\rm Z/H$')
plt.ylabel(r'$\rm a(\mu m)$')
plt.title(r'$\rm 100AU$')
plt.tick_params(pad=10,direction ='out')
plt.tick_params(length=6, width=1, which='major')
plt.tick_params(length=4, width=1, which='minor')

ax = plt.gca()
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="8%", pad=0.2)
cbar = plt.colorbar(myplt, cax=cax)
#plt.savefig('01m100au.eps',format='eps',dpi=300,bbox_inches='tight')
plt.show()
