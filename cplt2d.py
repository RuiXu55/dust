import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from pylab import *
from matplotlib import rc, font_manager
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable


au = '50'
s1 = '01'
out_ind   = [s1 +'m'+au+'auv_2.tab']
lab =['-','--',':']
############################################################################

name = 'dust_size_distribution'
f = open(out_ind[0],'r')
mark = True
while(mark):
 line = f.readline()
 data = line.split()
 if data[0]==name:
   mark = False
   numx = int(data[1])
   size = np.linspace(0,0,numx)
   for i in range(0,numx):
     line = f.readline()
     data = line.split()
     size[i] = float(data[0])
   break
f.close()

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
    if s==0:
      f = open(out_ind[l],'r')
    else:
      f = open(out_ind1[l],'r')

    numx = numx -1
    name = 'final'
    mark = True
    while(mark):
     line = f.readline()
     data = line.split()
     if data[0]==name:
       mark = False
       col = []
       numy = int(data[2])
       print ('numy',numy)
       print ('numx',numx)
       sum2 = np.ones((numx,numy-2))
       for i in range(0,numx):
         fi = f.readline()
         data = fi.split()
         col.append([])
         for k in range(0,numy):
           col[i].append(float(data[k]))
           if k>2 and k<numy-1:
             sum2[i,k-2] = abs(float(data[k])/(size[k-1]-size[k-3])*size[k-2]/neu[i])/2.0
           elif k==2:
             sum2[i,k-2] = abs(float(data[k])/(size[k-1]-size[k-2])*size[k-2]/neu[i])/2.0
           elif k==numy-1:
              sum2[i,k-2] = abs(float(data[k])/(size[k-2]-size[k-3])*size[k-2]/neu[i])/2.0
           print i,k-2,sum2[i,k-2]

######################## Plot the results #############################
Z = sum2.transpose()
extents = (0.0,px[numx-1],0.1,1e4)
cmap = cm.get_cmap('Greys_r',20)    # 11 discrete colors
cmap = cm.get_cmap('jet')    # 11 discrete colors
myplt = plt.imshow(Z,extent= extents,cmap=cmap,interpolation="bicubic",
      aspect='auto',origin='lower',norm=LogNorm(vmin=0.001, vmax=1.0))

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
plt.title(r'$\rm 50AU$')
plt.tick_params(pad=10,direction ='out')
plt.tick_params(length=6, width=1, which='major')
plt.tick_params(length=4, width=1, which='minor')

ax = plt.gca()
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="8%", pad=0.2)
cbar = plt.colorbar(myplt, cax=cax)
plt.savefig('01m50au.eps',format='eps',dpi=300,bbox_inches='tight')
plt.show()
