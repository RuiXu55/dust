import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
num = 200
xi = np.log(1000.0)/2.0
z = np.linspace(0,5,num)
alpha = np.linspace(0,0,num)
T = np.linspace(0,0,num)
for i in range(0,num):
  if z[i]<2.0:
    alpha[i] = 1e-4
  elif z[i]<4.0:
    alpha[i] = 1e-4*np.exp((z[i]-2.0)*xi)
  else:
    alpha[i] = 1e-4*np.exp(2.0*xi)

fig, ax1 = plt.subplots()
ax1.semilogy(z,alpha, 'k-',linewidth=2)
ax1.set_xlabel('Z/H')
# Make the y-axis label and tick labels match the line color.
ax1.set_ylabel(r'$\alpha_z$', color='k')
ax1.axis([0,5,1e-5,1e0])
for tl in ax1.get_yticklabels():
  tl.set_color('k')

T0 = 150.0*50**(-0.5)
for i in range(0,num):
  if z[i]<2.0:
    T[i] =T0
  elif z[i]<3.0:
    T[i] = T0 + (z[i]-2.0)*2.*T0
  else:
    T[i] = 3.*T0
ax2 = ax1.twinx()
ax2.plot(z,T, 'b-',linewidth=2)
ax2.axis([0,5,15,75])
ax2.set_ylabel('Temperature(K)', color='b')
for tl in ax2.get_yticklabels():
  tl.set_color('b')
'''
font = {'family' : 'sans-serif',
        'sans-serif' : 'Helvetica',
        'weight' : 'light',
        'size'   : 24}
matplotlib.rc('font', **font)
'''
matplotlib.rcParams.update({'font.size': 22})
plt.savefig('alpha.eps',format='eps',dpi=300,bbox_inches='tight')
plt.show()


