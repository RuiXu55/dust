import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


#out_ind ='1s40auc4'
out_ind ='ndif'
name = 'init_neutral_mass'
# find the name of the plotting
f = open(out_ind+'.tab','r')
while(1>0):
 line = f.readline()
 data = line.split()
 if data[0]==name:
   numx = int(data[1])
   numy = int(data[2])
   break
######################
neu = []
for i in range(0,numx):
  fi = f.readline()
  data = fi.split()
  neu.append(float(data[1]))
f.close()


name  = 'final'
# find the name of the plotting
f = open(out_ind+'.tab','r')
while(1>0):
 line = f.readline()
 data = line.split()
 if data[0]==name:
   numx = int(data[1])
   numy = int(data[2])
   break
######################

col = []
sum1 = np.linspace(0,0,numx) 
for i in range(0,numx):
  fi = f.readline()
  data = fi.split()
  col.append([])
  for k in range(0,numy):
    col[i].append(float(data[k]))
    if k>1:
      sum1[i] += float(data[k])

  ###### Plot results ##############
col = np.matrix(col)
col = col.T
px = np.linspace(0,0,numx)
py = np.linspace(0,0,numx)
for i in range(0,numx):
  px[i] = col[0,i]
  py[i] = col[1,i]/neu[i]

plt.plot(px,py,'k-',linewidth=2.0,label=r'$\alpha = \rm 10^{-4}$')
#plt.plot(px,py,'k-',linewidth=2.0,label=r'$1\mu m$')
f.close()

##################################################################
##################################################################
out_ind ='1s40auc3'
name = 'init_neutral_mass'
# find the name of the plotting
f = open(out_ind+'.tab','r')
while(1>0):
 line = f.readline()
 data = line.split()
 if data[0]==name:
   numx = int(data[1])
   numy = int(data[2])
   break
######################
neu = []
for i in range(0,numx):
  fi = f.readline()
  data = fi.split()
  neu.append(float(data[1]))
f.close()


name  = 'final'
# find the name of the plotting
f = open(out_ind+'.tab','r')
while(1>0):
 line = f.readline()
 data = line.split()
 if data[0]==name:
   numx = int(data[1])
   numy = int(data[2])
   break
######################

col = []
sum1 = np.linspace(0,0,numx) 
for i in range(0,numx):
  fi = f.readline()
  data = fi.split()
  col.append([])
  for k in range(0,numy):
    col[i].append(float(data[k]))
    if k>1:
      sum1[i] += float(data[k])

  ###### Plot results ##############
col = np.matrix(col)
col = col.T
px = np.linspace(0,0,numx)
py = np.linspace(0,0,numx)
for i in range(0,numx):
  px[i] = col[0,i]
  py[i] = col[1,i]/neu[i]

plt.plot(px,py,'k--',linewidth=2.0,label=r'$\alpha = \rm 10^{-3}$')
#plt.plot(px,py,'k--',linewidth=2.0,label=r'$1\mu m-1cm$')
f.close()

#################################################################
#################################################################
out_ind ='1s40auc5'
name = 'init_neutral_mass'
# find the name of the plotting
f = open(out_ind+'.tab','r')
while(1>0):
 line = f.readline()
 data = line.split()
 if data[0]==name:
   numx = int(data[1])
   numy = int(data[2])
   break
######################
neu = []
for i in range(0,numx):
  fi = f.readline()
  data = fi.split()
  neu.append(float(data[1]))
f.close()


name  = 'final'
# find the name of the plotting
f = open(out_ind+'.tab','r')
while(1>0):
 line = f.readline()
 data = line.split()
 if data[0]==name:
   numx = int(data[1])
   numy = int(data[2])
   break
######################

col = []
sum1 = np.linspace(0,0,numx) 
for i in range(0,numx):
  fi = f.readline()
  data = fi.split()
  col.append([])
  for k in range(0,numy):
    col[i].append(float(data[k]))
    if k>1:
      sum1[i] += float(data[k])

  ###### Plot results ##############
col = np.matrix(col)
col = col.T
px = np.linspace(0,0,numx)
py = np.linspace(0,0,numx)
for i in range(0,numx):
  px[i] = col[0,i]
  py[i] = col[1,i]/neu[i]
f.close()

#plt.plot(px,py,'k:',linewidth=2.0,label=r'$\rm 0.1\mu m-1cm $')
plt.plot(px,py,'k:',linewidth=2.0,label=r'$\rm \alpha = \rm 10^{-5} $')

#################################################################
out_ind ='1s40auv'
name = 'init_neutral_mass'
# find the name of the plotting
f = open(out_ind+'.tab','r')
while(1>0):
 line = f.readline()
 data = line.split()
 if data[0]==name:
   numx = int(data[1])
   numy = int(data[2])
   break
######################
neu = []
for i in range(0,numx):
  fi = f.readline()
  data = fi.split()
  neu.append(float(data[1]))
f.close()


name  = 'final'
# find the name of the plotting
f = open(out_ind+'.tab','r')
while(1>0):
 line = f.readline()
 data = line.split()
 if data[0]==name:
   numx = int(data[1])
   numy = int(data[2])
   break
######################

col = []
sum1 = np.linspace(0,0,numx) 
for i in range(0,numx):
  fi = f.readline()
  data = fi.split()
  col.append([])
  for k in range(0,numy):
    col[i].append(float(data[k]))
    if k>1:
      sum1[i] += float(data[k])

  ###### Plot results ##############
col = np.matrix(col)
col = col.T
px = np.linspace(0,0,numx)
py = np.linspace(0,0,numx)
for i in range(0,numx):
  px[i] = col[0,i]
  py[i] = col[1,i]/neu[i]
f.close()

plt.plot(px,py,'k-.',linewidth=2.0,label=r'$\rm \alpha \neq const $')


font = {'family' : 'sans-serif',
        'sans-serif' : 'Helvetica',
        'weight' : 'light',
        'size'   : 24}
matplotlib.rc('font', **font)

plt.axis([0,5,0,1])
plt.xlabel(r'$\rm Z/H $')
plt.ylabel(r'$\rm n_f/n_i$')
plt.title(r'$ \rm 40AU $')
plt.legend(loc=2)
fig = plt.gcf()
#plt.savefig(out_ind+'.eps',format='eps',dpi=300,bbox_inches='tight')
plt.savefig('40AUTc.eps',format='eps',dpi=300,bbox_inches='tight')
plt.show()
