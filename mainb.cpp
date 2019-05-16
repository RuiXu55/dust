/* This code calculates volatile molecules
 * evolution 
 */
 
#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <algorithm>
#include "header.hpp"

#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
using namespace boost::numeric::ublas;
using namespace std;

// template<class T> bool MInvBoost(const matrix<T>& InputMatrix, matrix<T>& InverseMatrix)
// Simplied header for matrices of doubles only, no other types allowed
void MInvBoost(matrix<double> InputMatrix, matrix<double> &InverseMatrix)
{
 // Create a duplicate of the input matrix
 int N = InputMatrix.size1();
 matrix<double> A = InputMatrix;
 // Create the permutation matrix
 typedef permutation_matrix<std::size_t> pmatrix;
 pmatrix P(N);
 // Assign the identity matrix to the inverse
 InverseMatrix.assign(identity_matrix<double>(N));
 //  LU factorization and substitution
 lu_factorize(A,P);
 lu_substitute(A,P,InverseMatrix);
 //  return true;
}

// dust size distribution
std::vector<double> Ds(double mins,double maxs,double 
  sigmad,double rho,int num,double xi){
  std::vector<double> Sd;
  double da,Dm;
  double a = 0.5;
  double sum = 0.0;
  for(int i=0;i<num;i++){
    if (num==1){
      da = mins;
    }else{
     da = mins*exp(i*log(maxs/mins)/(num-1.0));
    }
    Sd.push_back(pow(da,a)); 
    sum += Sd[i];
  }
  for(int j=0;j<num;j++){
    Sd[j] *= sigmad/sum;
  }
  return Sd;
}

// MAIN CODE STARTS HERE//
int main()
{
int Num,NDsize;
double Radius,Zmax,dZ,Evo_time,CFL,dt,tneu,tman;
double Sigma0,T0,Tm,Sigs,Tes,Sigma,Stick,dif;
double Dsize,Dmr,Drho,mins,maxs,xi,maxman,Mdust;
double Gsize,Gnr,Gm,Gmr,Gr,Nu0,Ed,ngr;
double sumrho,srho,rat,iso,sum,msiteden,sumH;
std::vector <double> Z,N_h,Ts,Vz,Sc,g1,
       n_neu,n_neu0,SigmaD,rhoz,T_hop,nlayer;
std::vector<Dust> Grain;
Dust ds;
myD dat,invd; 
ofstream f;

// Input parameter
f.open("1s100auv_omega.tab");
Radius   = 100.0;      // AU
mins     = 1.0;     // micron
NDsize   = 1;
Evo_time = 1.e6;       // Year
//Evo_time = 1.e4;       // Year
dif      = 1.0;

iso      = -1.0;
Num      = 100; 
Zmax     = 5.05;

Sigma0   = 500.0;
T0       = 150.0;
Sigs     = -1.0;
Tes      = -0.5;
Dmr      = 1e-2;
maxs     = 1e4;     // micron
Drho     = 2.0;
Gsize    = 1e-5;
Gnr      = 1e-4;
Gm       = 28;
Ed       = 1150.0;

CFL      = 0.002;
dt       = 1e10;    // initial
Stick    = 1.0;

dZ       = Zmax/(Num-1.0);
Tm       = T0*pow(Radius,Tes);
Sigma    = Sigma0*pow(Radius,Sigs);
Mdust    = Sigma*Dmr;
Nu0      = sqrt(2.50903e22*Ed/Gm);
double sum0 = 0.0;
double time = 0.0;

//cout<<"Simulation Starts!"<<endl;

// vertical density profile
for(int i=0;i<Num;i++){
  rhoz.push_back(0.0);
  Z.push_back(i*dZ);
}

if(iso>0.0){
  for (int i=0;i<Num;i++){
  rhoz[i] = Hrho(Radius,Tm,Z[i],Sigma);
  }
}else{
  srho = 0.0;
  for (int j=0;j<Num;j++){
    sumrho =0.0;
    for (int i=0;i<=j;i++){
      if(i==0){
      sumrho -= Z[i]*dZ/pow(Cs(Z[i]),2);
      }else{
      sumrho -= (Z[i]*dZ+(pow(Cs(Z[i]),2)
         -pow(Cs(Z[i-1]),2)))/pow(Cs(Z[i]),2);
      }
    }
  rhoz[j] = exp(sumrho);
  srho += rhoz[j]*dZ; 
  }
  rat = Sigma/2.0/H(Radius,Tm)/srho;
  double sumsigma;
  sumsigma = 0.0;
  for(int i=0;i<Num;i++){
    rhoz[i] = rhoz[i]*rat;
    sumsigma += rhoz[i]*dZ;
    //cout<<" rhoz "<<rhoz[i]<<endl;
  }
  //cout<<"sigma="<<sumsigma<<endl;
}


/* GAS PROPERTY */
for (int i=0;i<Num;i++){
  // hydrogen # density
  N_h.push_back(Hn(rhoz[i]));
  // gas density
  Gr = Grho(rhoz[i],Gnr,Gm);
  // non-dimensionless time
  Ts.push_back(Tausg(Radius,Tm,Z[i],rhoz[i],Gr,Gsize));
  // Schmit #
  Sc.push_back(1.+pow(Ts[i],2));
  //Vz.push_back(-min(Cs(Radius,Tm,Z[i]),Ts[i]*Z[i]));
  Vz.push_back(-Ts[i]*Z[i]/Sc[i]);
  // gas # density
  n_neu.push_back(Gn(rhoz[i], Gnr));
  n_neu0.push_back(Gn(rhoz[i], Gnr));
  // total gas
  if(i>0 && i<Num-1){
    sum0 += n_neu[i];
  }
}

// reflection boundary condition
Vz[0]    = -Vz[1];
Vz[Num-1]= -Vz[Num-2];
// diffusion coefficient
for(int i=0;i<Num-1;i++){
  g1.push_back((N_h[i+1]*Al(Z[i+1])*Cs(Z[i+1])/Sc[i+1]+
          N_h[i]*Al(Z[i])*Cs(Z[i])/Sc[i])/2.0);
}
cout<<"Gas initialized!"<<endl;

// dust properties
// dust mass at different size
SigmaD= Ds(mins,maxs,Dmr,Drho,NDsize,xi);

for(int j=0;j<NDsize;j++){
  ds = Dust();
  if (NDsize==1){
    ds.size = mins;
  }else{
   ds.size = mins*exp(j*log(maxs/mins)/(NDsize-1));
  }

  ds.mr=SigmaD[j];
  for(int i=0;i<Num;i++){
    ds.Ts.push_back(Taus(Radius,Tm,Z[i],rhoz[i],Drho,ds.size));
    ds.Sc.push_back(1.+ pow(ds.Ts[i],2));
    ds.Vz.push_back(-ds.Ts[i]*Z[i]/ds.Sc[i]);
    ds.man.push_back(0.0);

    /* evolution timestep */
    if(i>0){
      dt = min(dt,CFL*dZ/max(ds.Vz[i],-ds.Vz[i]));
      dt = min(dt,0.1*pow(dZ,2)/(Al(Z[i])*Cs(Z[i])/ds.Sc[i]));
    }
  }
  // reflection boundary condition
  ds.Vz[0] = -ds.Vz[1]; 
  ds.Vz[Num-1]= -ds.Vz[Num-2];
  // diffusion coefficient
  for(int i=0;i<Num-1;i++){
   ds.g1.push_back((N_h[i+1]*Al(Z[i+1])*Cs(Z[i+1])/ds.Sc[i+1]+
      N_h[i]*Al(Z[i])*Cs(Z[i+1])/ds.Sc[i])/2.0);
  }

  // dust Steady state vertical distribution.
  double mratio = 0.0;
  double Thydro = 0.0;
  for ( int i=0;i<Num;i++){
     ds.fd.push_back(integral(Radius,Tm,dZ,i,rhoz[i],Drho,ds.size));
     mratio += ds.fd[i]*N_h[i]*dZ;
     Thydro += N_h[i]*dZ;
  }
  for (int i=0;i<Num;i++){
    ds.fd[i] *= ds.mr/(mratio/Thydro);
  }
  
  // adsorp/desorp rate
  for (int i=0;i<Num;i++){
    ds.ad.push_back(2.844e-20*Stick/ds.size*
     sqrt(Td(Tm,Z[i])/50.0)/sqrt(Gm/28.0)
     *(ds.fd[i]/1.e-4)*
     N_h[i]*dZ/Omega(Radius));

    ds.des.push_back(Nu0*exp(-Ed/Td(Tm,Z[i]))
     /Omega(Radius));

  }

  /* invert matrix of chemical reactions */
  for (int i=0;i<Num;i++){
    dat.d1 = 1.0 + dt*ds.ad[i];
    dat.d2 = -ds.ad[i]*dt;
    dat.d3 = -ds.des[i]*dt;
    dat.d4 = 1.0+ dt*ds.des[i];
    double coef = 1.0/(dat.d1*dat.d4-dat.d2*dat.d3);
    invd.d1 = coef*dat.d4;
    invd.d2 = -coef*dat.d2;
    invd.d3 = -coef*dat.d3;
    invd.d4 = coef*dat.d1;
    ds.inv.push_back(invd);
    }
  Grain.push_back(ds);
}
cout<<"Dust Initialize Successful!"<<endl;

/* Output initial value */
// output grain mass size distribution
f<<"dust_size_distribution "<<NDsize<<" "<<2<<endl;
for(int i=0;i<NDsize;i++){
  f<<setw(15)<<Grain[i].size<<setw(15)<<Grain[i].mr<<endl;;
}
// initial gas number density
f<<"init_neutral_mass "<<Num<<" "<<2<<endl;
for(int i=0;i<Num;i++){
  f<<setw(15)<<Z[i]<<setw(15)<<n_neu[i]<<endl;;
}
// initial dust vertical profile
f<<"dust_vert_distribution "<<Num<<" "<<NDsize+1<<endl;
for(int i=0;i<Num;i++){
  f<<setw(15)<<Z[i];
  for (int j=0;j<NDsize;j++){
  f<<setw(15)<<Grain[j].fd[i];
  }
  f<<endl;
}

/*
 * Time Evolution starts
 */

std::vector<double> adflux,difflux;
for(int i=0;i<Num-1;i++){
  adflux.push_back(0.0);
  difflux.push_back(0.0);
}
for(int i=0;i<NDsize;i++){
  nlayer.push_back(0.0);
}


double outtime = 0.0;
f<<"Evolution starts!"<<endl;
while (time<Evo_time){
  time += dt/Omega(Radius)/3.154e7;
  outtime += dt/Omega(Radius)/3.154e7;
  //time += dt;
  //outtime += dt;
  if (outtime>1000.0){
    sum = 0.0;
    for (int j=1;j<Num-1;j++){
     sum += n_neu[j];
     for (int i=0;i<NDsize;i++){
        sum += Grain[i].man[j];
      }
    }
    cout<<"Time="<<time<<" Sum= "<<sum<<endl;
    f<<"timestep "<<Num<< "  "<<NDsize+2<<endl;
    for(int i=0;i<Num;i++){
      f<<setw(15)<<Z[i]<<setw(15)<<n_neu[i];
      for(int j=0;j<NDsize;j++){
        f<<setw(15)<<Grain[j].man[i];
      }
      f<<endl;
      }
    outtime = 0.0;
    }
   
  
  // updating desorption and adsorption
  for (int i=0;i<NDsize;i++){
  for (int j=0;j<Num;j++){
    // mantle species is single layer or multiple layers 
    msiteden = 1.5e15*4.0*3.1416*pow(Grain[i].size/1.e4,2);
    ngr = N_h[j]*Grain[i].fd[j]*1.6726e-24*2.34/
      (4.0/3.0*3.1416*Drho*pow(Grain[i].size/1.e4,3));
    maxman = msiteden*ngr;
    nlayer[i] = Grain[i].man[j]/maxman;

    tman = Grain[i].man[j]+dt*(n_neu[j]*Grain[i].ad[j]-
      Grain[i].man[j]*Grain[i].des[j]/max(1.0,nlayer[i]));
    //tman = Grain[i].man[j]+n_neu[j]*Grain[i].ad[j]-
    //  Grain[i].man[j]*Grain[i].des[j];
    if(tman<0.0)
      tman = 0.0;
    if(tman>n_neu[j]+Grain[i].man[j])
      tman = n_neu[j] + Grain[i].man[j];
    n_neu[j] += Grain[i].man[j]-tman;
    Grain[i].man[j] = tman;
    }

    // boundary condition 
    Grain[i].man[Num-1] = Grain[i].man[Num-2];
    Grain[i].man[0] = Grain[i].man[1];
  }
  // boundary condition
  n_neu[Num-1] = n_neu[Num-2];
  n_neu[0] = n_neu[1];
  

  // Convection and diffusion
  if (dif>0.0){ 
       for(int j=1;j<Num-1;j++){
         adflux[j]  = 0.0*n_neu[j+1]*Vz[j+1]; 
         difflux[j] = -g1[j]*(n_neu[j+1]/N_h[j+1]-n_neu[j]/N_h[j])/dZ;
      }
      adflux[0] = adflux[Num-2] = 0.0;
      difflux[0] = difflux[Num-2]= 0.0;
      for(int j=1;j<Num-1;j++){
        n_neu[j] += dt*(-adflux[j]-difflux[j]
                    +adflux[j-1] +difflux[j-1])/dZ;
      }
      // boundary condition
      n_neu[Num-1] = n_neu[Num-2];
      n_neu[0] = n_neu[1];

    // n_man should follow the same adv/diffu equation as dust
    for (int i=0;i<NDsize;i++){
        for(int j=1;j<Num-1;j++){
           adflux[j] = Grain[i].man[j+1]*Grain[i].Vz[j+1]; 
           difflux[j]= -Grain[i].g1[j]*(Grain[i].man[j+1]/N_h[j+1]-
                       Grain[i].man[j]/N_h[j])/dZ;
         }
         adflux[0] = adflux[Num-2]= 0.0;
         difflux[0] = difflux[Num-2] = 0.0;

        for(int j=1;j<Num-1;j++){
          Grain[i].man[j] += dt*(-adflux[j]-difflux[j]+
                             adflux[j-1]+difflux[j-1])/dZ;
        }
          Grain[i].man[Num-1] = Grain[i].man[Num-2];
          Grain[i].man[0] = Grain[i].man[1];
      }
  }

}

f<<"final "<<Num<< "  "<<NDsize+2<<endl;
for(int i=0;i<Num;i++){
  f<<setw(15)<<Z[i]<<setw(15)<<n_neu[i];
  for(int j=0;j<NDsize;j++){
    f<<setw(15)<<Grain[j].man[i];
    }
  f<<endl;
}
cout<<"Evolution succeed!"<<endl;

return 0;
}
