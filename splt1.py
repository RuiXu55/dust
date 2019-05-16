import sys
import numpy as np
import matplotlib
import matplotlib.font_manager
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'

font = {'family' : 'sans-serif',
        'sans-serif' : 'Helvetica',
        'weight' : 'light',
        'size'   : 18}
matplotlib.rc('font', **font)


out_ind  = '1s50auv_2e6.tab'
out_ind1   = '1s50auv_2e6.tab'

Iter = 1000
wid = 2.0

############################################################################
name = 'init_neutral_mass'
# find the name of the plotting
f = open(out_ind,'r')
mark = True
while(mark):
 line = f.readline()
 data = line.split()
 if data[0]==name:
   mark = False
   numx = int(data[1])
   px = np.linspace(0,0,numx)
   py = np.linspace(0,0,numx)
   py0 = np.linspace(0,0,numx)
   neu = np.linspace(0,0,numx)
   for i in range(0,numx):
     line = f.readline()
     data = line.split()
     px[i] = float(data[0])
     neu[i] = float(data[1])
     py[i] = 1.0
     py0[i] = 0.0
   break

f1 = open(out_ind1,'r')
while(1>0):
 line1 = f1.readline()
 data1 = line1.split()
 if data1[0]==name:
   numx1 = int(data1[1])
   px1 = np.linspace(0,0,numx1)
   py1 = np.linspace(0,0,numx1)
   py2 = np.linspace(0,0,numx1)
   neu1 = np.linspace(0,0,numx1)
   for i in range(0,numx1):
     line = f1.readline()
     data = line.split()
     px1[i] = float(data[0])
     neu1[i] = float(data[1])
     py1[i] = 1.0
     py2[i] = 0.0
   break
print ('numx1=',numx1)

name = 'dust_vert_distribution'
f = open(out_ind,'r')
mark = True
while(mark):
 line = f.readline()
 data = line.split()
 if data[0]==name:
   mark = False
   numy = int(data[2])
   break
f.close()

f = open(out_ind1,'r')
mark = True
while(mark):
 line = f.readline()
 data = line.split()
 if data[0]==name:
   mark = False
   numy1 = int(data[2])
   break
f.close()
f1.close()


#####################################
name = 'timestep'
f = open(out_ind,'r')
f1 = open(out_ind1,'r')
while(1>0):
 line = f.readline()
 data = line.split()
 if data[0]==name:
   break
while(1>0):
 line = f1.readline()
 data = line.split()
 if data[0]==name:
   break

for n in range(0,Iter):
  print ('n=',n)
  sum0 = np.linspace(0,0,numx) 
  sum1 = np.linspace(0,0,numx1) 
  ############################
  for i in range(0,numx):
    fi = f.readline()
    data = fi.split()
    px[i] = float(data[0])
    py[i] = float(data[1])/neu[i]
    for j in range(2,numy+1):
      sum0[i] += float(data[j])
    py0[i]= sum0[i]/neu[i]
    #print px[i],py[i],py0[i]
  #############################
  for i in range(0,numx1):
    fi1 = f1.readline()
    data1 = fi1.split()
    px1[i] = float(data1[0])
    py1[i] = float(data1[1])/neu1[i]
    for j in range(2,numy1+1):
      sum1[i] += float(data1[j])
    py2[i]= sum1[i]/neu1[i]

  if n<Iter-1:
    fi = f.readline()
    fi1 = f1.readline()
  ####################################################
  ind1 = 0
  if n==9 or n==99 or n==999:
    ind1 += 1
    if n==9:
      plt.plot(px,py,'k-',linewidth=wid,label=r'$\rm gas$')
      plt.plot(px,py0,'b-',linewidth=wid,label=r'$\rm mantle$')
      #plt.semilogy(px,py,'k-',linewidth=wid,label=r'$\rm gas$')
      #plt.semilogy(px,py0,'b-',linewidth=wid,label=r'$\rm mantle$')
    else:
      #plt.semilogy(px,py,'k-',linewidth=wid)
      #plt.semilogy(px,py0,'b-',linewidth=wid)
      plt.plot(px,py,'k-',linewidth=wid)
      plt.plot(px,py0,'b-',linewidth=wid)
    plt.plot(px1,py1,'k--',linewidth=wid)
    plt.plot(px1,py2,'b--',linewidth=wid)
    #plt.annotate(str(ind1), xy=(px[20], py[20]), 
    #                xytext=(20, 10), textcoords='offset points', va='center',
    #                            arrowprops=dict(arrowstyle='->'))
    font = {'family' : 'sans-serif',
            'sans-serif' : 'Helvetica',
            'weight' : 'light',
            'size'   : 20}
    matplotlib.rc('font', **font)
    plt.axis([0,5,0.0,1.6])
    plt.xlabel(r'$\rm Z/H$')
    plt.ylabel(r'$\rm n/n_{initial}$')
    #plt.title(r'$\rm 100AU$')
    plt.legend(loc=0)

'''
#plt.text(1.7,1.05,r'1',fontsize=20)
plt.text(1.75,0.85,r'1',fontsize=20)
plt.text(1.75,0.5,r'2',fontsize=20)
plt.text(1.75,0.07,r'3',fontsize=20)
#plt.text(2.3,1.05,r'1',fontsize=20)
plt.text(2.4,0.85,r'1',fontsize=20)
plt.text(2.4,0.35,r'2',fontsize=20)
plt.text(2.4,0.07,r'3',fontsize=20)
'''
f.close()
f1.close()
plt.savefig('evo.pdf',format='pdf',dpi=500,bbox_inches='tight')
plt.show()
