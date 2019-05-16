#include <iostream>
#include <math.h>

// Omega
double Omega(double r){
  //return sqrt(1.3275e26/pow(r*1.496e13,3));
  return 1.99126e-7*pow(r,-3.0/2.0);
}

// Disk Temperature
double Td(double Tm, double z){
  //return Tm;
  if (z<2.0)
    return Tm;
  else if (z<3.0)
    return Tm+2.0*Tm*(z-2.0);
  else
    return 3.0*Tm;
}

// midplane scale height (cs/omega)
double H(double r, double Tm){
  return sqrt(3.5275e7*Tm)/Omega(r);
}

// sound speed ( normalized)
double Cs(double z){
  if (z<2.0)
    return 1.0;
  else if (z<3.0)
    return sqrt(1.0+2.0*(z-2.0));
  else
    return sqrt(3.0);
}

// gas thermal speed (cgs)
double Cs1(double Tm,double z){
  return sqrt(3.528e7*Td(Tm,z));
}


// co thermal speed (cgs)
double Vth(double r,double Tm,double z){
  return sqrt(7.507e6*Td(Tm,z));
}

// gas thermal speed (cgs)
double Vthg(double r,double Tm,double z){
  return sqrt(8.983e7*Td(Tm,z));
}

// isothermal hydrogen density
double Hrho(double r, double Tm,double z,double sigma){
  return sigma/(2.50663*H(r,Tm))*exp(-pow(z,2)/2.0);
}

// hydrogen # density
double Hn(double rhoz){
  return rhoz/1.6726e-24;
}

// gas # density
double Gn(double rhoz, double Gnr){
  return Gnr*Hn(rhoz);
}

// gas mass density
double Grho(double rhoz,double Gnr,double Gm){
  return rhoz*Gnr*Gm;
}

// Dz
double Al(double z){
  //return 1.e-2;
  double xi = log(1000.0)/2.0;
  if(z<2.0)
    return 1e-4;
  else if (z<4.0)
    return 1e-4*exp((z-2.0)*xi);
  else
    return 1e-4*exp(2.0*xi);
}

// non-dimension stopping time
double Taus(double r, double Tm,double z,double rhoz,double Drho,double Gsize){  
 return Drho/rhoz*(Gsize*1.e-4)*Omega(r)/Cs1(Tm,z);
 //return Drho/rhoz*(Gsize*1e-4)*Omega(r)/Vthg(r,Tm,z);
}

double Tausg(double r, double Tm,double z,double rhoz,double Drho,double Gsize){ 
 //return Drho/rhoz*(Gsize*1e-4)*Omega(r)/Vthg(r,Tm,z);
 return Drho/rhoz*(Gsize*1.e-4)*Omega(r)/Cs1(Tm,z);
}

// Schimt number
double Sc(double r, double Tm, double z, double rhoz,double Drho,double size){ 
 return 1.+pow(Taus(r,Tm,z,rhoz,Drho,size),2);
}

double integral(double r, double Tm, double dz,int ind, double rhoz,double rho, double size)
{
     double sum,z,Ts,Sc;
     sum = 0.0;
     for( int i=0;i<ind+1;i++){
       //z = (i+0.5)*dz;
       z  = i*dz;
       Ts = Taus(r,Tm,z,rhoz,rho,size);
       Sc = 1.0 + pow(Ts,2); 
       //sum += dz*(-Sc*Ts*z/Al(z));         
       sum += dz*(-Ts*z/Al(z)/Cs(z));         
     } 
      return exp(sum);
}


