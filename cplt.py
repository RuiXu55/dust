import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

au = '50'
s1 = '1'
out_ind   = ['1s'+au+'auv_2e6.tab',
             '1s'+au+'auc4.tab',
             '1s'+au+'auc2.tab']
#out_ind   = ['1m'+au+'auv.tab',
#             '01m'+au+'auv.tab']
lab =['-','--',':']
#labels=[r'$\alpha_z \neq \rm const$',r'$\alpha_z = 10^{-4}$',r'$\alpha_z = 10^{-2}$']
#labels=[r'$a=10^{-1}-10^4 \mu m$',r'$a=1.0-10^4\mu m$']
wid = 2.0
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
    if s==0:
      f = open(out_ind[l],'r')
    else:
      f = open(out_ind1[l],'r')

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
       sum1 = np.linspace(0,0,numx) 
       for i in range(0,numx):
         fi = f.readline()
         data = fi.split()
         col.append([])
         for k in range(0,numy):
           col[i].append(float(data[k]))
           if k>1:
              sum1[i] += float(data[k])
       col = np.matrix(col)
       col = col.T
       py = np.linspace(0,0,numx)
       py0 = np.linspace(0,0,numx)
       for i in range(0,numx):
          px[i] = col[0,i]
          py[i] = col[1,i]/neu[i]
          py0[i] = sum1[i]/neu[i]
    py[0] = 2.*py[1]-py[2]
    py0[0] = 2.*py0[1]-py0[2]
    ######################
    if l==0:
      plt.plot(px,py,'k-',linewidth=wid,label=r'$\rm gas $')
      plt.plot(px,py0,'b-',linewidth=wid,label=r'$\rm mantle$')
    else:
      plt.plot(px,py,'k',linestyle=lab[l],linewidth=wid)
      plt.plot(px,py0,'b',linestyle=lab[l],linewidth=wid)
    '''
    plt.plot(px,py,'k',linestyle=lab[l],linewidth=wid,label=labels[l])
    plt.plot(px,py0,'b',linestyle=lab[l],linewidth=wid)
    '''

    font = {'family' : 'sans-serif',
            'sans-serif' : 'Helvetica',
            'weight' : 'light',
            'size'   : 24}
    matplotlib.rc('font', **font)
    plt.axis([0.,5,0.,1.6])
    plt.xlabel(r'$\rm Z/H$')
    plt.ylabel(r'$\rm n/n_{initial}$')
    plt.title(r'$\rm 100AU$')
    #plt.text(1.0,1.3,txt, fontsize=18)
    matplotlib.rc('font', **font)
    plt.legend(loc=1)
    f.close()
fig = plt.gcf()
fi = 's'+au+'AU.eps'
plt.savefig(fi,format='eps',dpi=300,bbox_inches='tight')
plt.show()

