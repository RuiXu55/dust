import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


out_ind ='tst'

name = 'init_neutral_mass'
# find the name of the plotting
f = open(out_ind+'.tab','r')
while(1>0):
 line = f.readline()
 data = line.split()
 if data[0]==name:
   numx = int(data[1])
   numy = int(data[2])
   print( 'numx',numx,'numy',numy)
   break
######################
neu = []
Z = []
for i in range(0,numx):
  fi = f.readline()
  data = fi.split()
  neu.append(float(data[1]))
  Z.append(float(data[0]))


name = 'final'
# find the name of the plotting
while(1>0):
 line = f.readline()
 data = line.split()
 if data[0]==name:
   numx = int(data[1])
   numy = int(data[2])
   print( 'numx',numx,'numy',numy)
   break
######################
man = []
gas = []
Z = []
for i in range(0,numx):
  fi = f.readline()
  data = fi.split()
  gas.append(float(data[1])/neu[i])
  man.append(float(data[2])/neu[i])
  Z.append(float(data[0]))
plt.plot(Z,gas,'k-')
plt.plot(Z,man,'r-')
plt.axis([0,5,0,1.5])
plt.show()
f.close()


