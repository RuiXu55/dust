import numpy as np
from scipy import integrate

x = np.log(1000.0)/2.0

Num  = 100
Z = np.linspace(0,5,Num)
dZ = Z[-1]/(Num-1)
Radius = 100
AU = 1.496e13

Sigma0 = 500
T0 = 150
Tdisk = T0*Radius**(-0.5)

rho = 2.0 
mr = 1e-2

Sigma    = Sigma0*Radius**(-1.0)
GM       = 6.67e-8*1.99e33
Omega    = np.sqrt(GM/(Radius*AU)**3)
H        = np.sqrt(3.5275*Tdisk)/Omega

def integrand(z,Rho1,Size1,ind):
   # tem
   if z<2.0:
    Tm = Tdisk
   elif z<3.0:
    Tm = Tdisk+2.0*Tdisk*(z-2.0)
   else:
    Tm = 3.0*Tdisk
   Tm = Tdisk
   Vth      = np.sqrt(7.507e6*Tm)
   print ('Vth',Vth)

   if z<2.0:
     al = 1e-4
   elif z<4.0:
     al  = 1e-4*np.exp((z-2.0)*x)
   else:
     al  = 1e-4*np.exp((4.0-2.0)*x)

   Rhoh = Sigma/(np.sqrt(2.0*np.pi)*H)*\
          np.exp(-z**2/2.0)
   Ts = Rho1/Rhoh*(Size1*1e-4)/Vth*Omega
   Sc = 1.+Ts**2
   return -Sc*Ts*z/al

def profile(rho,size,ra,ind):
 add_m   = 0
 m_ratio = []
 for i in range(0,Num):
   ans,err = integrate.quad(integrand,0,Z[i],args=(rho,size,ind))
   m_ratio =np.append(m_ratio,np.exp(ans))
   add_m += m_ratio[i]*N_h[i]*dZ

 const = ra/(add_m/np.sum(N_h*dZ))
 return m_ratio*const

if __name__ == "__main__":
  import matplotlib.pyplot as plt
  import matplotlib
  ind = np.array([0,1,2])
  size= np.array([0.1,1.0])
  colors = ['r','k','b']
  line = ['-','--',':']
  labels = [r'$\rm 0.1\mu m$',r'$\rm 10\mu m$']
  for i in range(0,3):
    for j in range(0,2):
      f_d = profile(rho,size[j],mr,ind[i])
      if i==0:
        plt.plot(Z,f_d*100.0,color=colors[j],linestyle=line[i],linewidth=2.0,label=labels[j])
      else:
        plt.plot(Z,f_d*100,color=colors[j],linestyle=line[i],linewidth=2.0)

  font = {'family' : 'sans-serif',
          'sans-serif' : 'Helvetica',
          'weight' : 'light',
          'size'   : 24}
  matplotlib.rc('font', **font)
  plt.legend(loc=0)
  plt.xlabel(r'$\rm Z/H$')
  plt.ylabel(r'$\rm \rho_d/\rho_g(\times 100)$')
  plt.title(r'$\rm 100 AU$') 
  #plt.savefig('100AU.eps',format='eps',dpi=300,bbox_inches='tight')
  plt.show()

