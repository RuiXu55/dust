import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
out_ind =['1s100auv','1s100auc4','1s100auc2']
out_ind1 =['01s100auv','01s100auc4','01s100auc2']
labels = [r'$\rm 1\mu m$',r'$\rm 0.1\mu m$']

name = 'dust_vert_distribution'
lab = ['-','-.','--']
col = ['k','r']
for x in range(0,3):
 print ('x=',x)
 for y in range(0,2):
   print ('y=',y)
   if y==0:
     f = open(out_ind[x]+'.tab','r')
   else:
     f = open(out_ind1[x]+'.tab','r')
   while(1>0):
     line = f.readline()
     data = line.split()
     if data[0]==name:
       numx = int(data[1])
       numy = int(data[2])
       print( 'numx',numx,'numy',numy)
       dust = np.linspace(0,0,numx)
       Z = np.linspace(0,0,numx)
       for i in range(0,numx):
         fi = f.readline()
         data = fi.split()
         Z[i] = float(data[0])
         dust[i] = float(data[1])
       if x==0: 
         plt.plot(Z,dust*100.0,linestyle=lab[x],
             color=col[y],linewidth=2.0,label=labels[y])
       else:
         plt.plot(Z,dust*100.0,linestyle=lab[x],
             color=col[y],linewidth=2.0)
       #plt.plot(Z,dust)
       break

matplotlib.rcParams.update({'font.size': 22})
plt.legend(loc=0)
plt.xlabel(r'$\rm Z/H$')
plt.ylabel(r'$\rm \rho_d/\rho_g(\times 100)$')
plt.title(r'$\rm 50 AU$') 
plt.axis([0,5,0,1.5])
#plt.savefig('dust.eps',format='eps',dpi=300,bbox_inches='tight')
plt.show()
f.close()


