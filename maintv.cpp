#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <algorithm>
#include "header.hpp"
using namespace std;

// dust size distribution
vector<double> Ds(double mins,double maxs,double sigmad,double rho,int num,double xi){
  vector<double> Sd;
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
double Dsize,Dmr,Drho,mins,maxs,xi;
double Gsize,Gnr,Gm,Gmr,Gr,Nu0,Ed;
vector <double> Z,N_h,T_hop,Ts,Vz,Sc,g1,n_neu,SigmaD;
vector<Dust> Grain;
myD dat,invd; 
ofstream f,f1;

// Input parameter
f.open("1s100auv_tv.tab");

dif      = 1.0;
Radius   = 100.0;      // AU
Num      = 120; 
Zmax     = 5.1;
Evo_time = 5e6;       // Year

Sigma0   = 500.0;
T0       = 150.0;
Sigs     = -1.0;
Tes      = -0.5;
Dmr      = 1e-2;
mins     = 1.0;     // micron
maxs     = 1e4;     // micron
NDsize   = 1;
Drho     = 2.0;
Gsize    = 1e-5;
Gnr      = 1e-4;
Gm       = 28;
Ed       = 1150.0;

CFL      = 0.4;
dt       = 1e10;    // initial
Stick    = 1.0;

dZ       = Zmax/(Num-1.0);
Tm       = T0*pow(Radius,Tes);
Sigma    = Sigma0*pow(Radius,Sigs);
Nu0      = sqrt(2.5118e22*Ed/Gm);
cout<<"Nu0"<<Nu0<<endl;
//return 0;

// height and  vertical density
vector<double> rhoz;
for(int i=0;i<Num;i++){
  Z.push_back(i*dZ);
}
//cout<<" Simulation Starts!"<<endl;
// vertical density
double sumrho,srho,rat;
for(int i=0;i<Num;i++) rhoz.push_back(0.0);
srho = 0.0;
for (int j=0;j<Num;j++){
  sumrho =0.0;
  for (int i=0;i<=j;i++){
    if(i==0){
    sumrho -= Z[i]*dZ/pow(Cs(Radius,Tm,Z[i]),2);
    }else{
    sumrho -= (Z[i]*dZ+pow(Cs(Radius,Tm,Z[i]),2)
        -pow(Cs(Radius,Tm,Z[i-1]),2))/pow(Cs(Radius,Tm,Z[i]),2);
    }
  }
rhoz[j] = exp(sumrho);
srho += rhoz[j]*dZ; 
}
rat = Sigma/2.0/H(Radius,Tm)/srho;
for(int i=0;i<Num;i++){
  rhoz[i] = rhoz[i]*rat;
  //rhoz[i] = Hrho(Radius,Tm,Z[i],Sigma);
  //cout<<rhoz[i]<<"  "<<Hrho(Radius,Tm,Z[i],Sigma)<<endl;
}

//cout<<" Simulation Starts!"<<endl;

// gas property
for (int i=0;i<Num;i++){
  // hydrogen # density
  N_h.push_back(Hn(rhoz[i]));
  // hopping time
  T_hop.push_back(exp(Ed/Td(Tm,Z[i]))/Nu0);
  // gas density
  Gr = Grho(rhoz[i],Gnr,Gm);
  // non-dimensionless time
  Ts.push_back(Taus(Radius,Tm,Z[i],rhoz[i],Gr,Gsize));
  // convection vel
  Vz.push_back(-min(Cs(Radius,Tm,Z[i]),Ts[i]*Z[i]));
  // Schmit #
  Sc.push_back(1.+pow(Ts[i],2));
  // gas # density
  n_neu.push_back(Gn(rhoz[i], Gnr));
}

// reflection boundary condition
Vz[0]    = -Vz[1];
Vz[Num-1]= -Vz[Num-2];
// diffusion coefficient
for(int i=0;i<Num-1;i++){
  g1.push_back((N_h[i+1]*Al(Z[i+1])/Sc[i+1]+
          N_h[i]*Al(Z[i])/Sc[i])/2.0);
}

// total gas
double sum0 = 0.0;
for (int j=1;j<Num-1;j++){
  sum0 += n_neu[j];
}
cout<<"Init Gas = "<<sum0<<endl;

//cout<<"Gas initialized!"<<endl;


// time step
Dust ds;
ds = Dust();
ds.size = mins;
for(int i=0;i<Num;i++){
  ds.Ts.push_back(Taus(Radius,Tm,Z[i],rhoz[i],Drho,ds.size));
  ds.Vz.push_back(-min(Cs(Radius,Tm,Z[i]),ds.Ts[i]*Z[i]));
  ds.Sc.push_back(1.+ pow(ds.Ts[i],2));
}
for(int i=1;i<Num;i++){
  dt = min(dt,CFL*dZ/max(ds.Vz[i],-ds.Vz[i]));
  dt = min(dt,0.1*pow(dZ,2)/(Al(Z[i])/ds.Sc[i]));
}
dt = dt/10.0;


// dust properties
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
    ds.Vz.push_back(-min(Cs(Radius,Tm,Z[i]),ds.Ts[i]*Z[i]));
    ds.Sc.push_back(1.+ pow(ds.Ts[i],2));
    ds.man.push_back(0.0);
  }

  // reflection boundary condition
  ds.Vz[0] = -ds.Vz[1]; 
  ds.Vz[Num-1]= -ds.Vz[Num-2];

  // diffusion coefficient
  for(int i=0;i<Num-1;i++){
   ds.g1.push_back((N_h[i+1]*Al(Z[i+1])/ds.Sc[i+1]+
      N_h[i]*Al(Z[i])/ds.Sc[i])/2.0);
  }
  // dust Steady state distribution.
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
    ds.ad.push_back(2.19e-21*Stick/ds.size*
     sqrt(Td(Tm,Z[i])/50.0)/sqrt(Gm/300.0)
     *(ds.fd[i]/1.e-4)*
     N_h[i]/Omega(Radius));

    ds.des.push_back(1.0/(T_hop[i]*Omega(Radius)));
    /*
    for (int x =0;x<1000;x++){
    double temx = 15.0*(1.0+x/1000.0);
    cout<<" tempx "<<temx<<endl;
    cout<<"ad "<<2.19e-21*Stick/ds.size*sqrt(temx/50.0)/sqrt(Gm/300.0)*(ds.fd[i]/1.e-4)*N_h[i]/Omega(Radius)<<endl;
    double hopp =exp(Ed/temx)/Nu0;
    cout<<"des "<<1.0/(hopp*Omega(Radius))<<endl;
    }
    return 0;
    */
  }

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

//cout<<" Dust Initialize Successful!"<<endl;
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

//cout<<"timestep="<<dt<<endl;

//#####################################################################
//######### Time Evolution starts
double sum;
double time = 0.0;
vector<double> adflux;
vector<double> difflux;
for(int i=0;i<Num-1;i++) adflux.push_back(0.0);
for(int i=0;i<Num-1;i++) difflux.push_back(0.0);


double outtime = 0.0;
f<<"evolution"<<endl;
while (time<Evo_time){
  time += dt/Omega(Radius)/3.154e7;
  outtime += dt/Omega(Radius)/3.154e7;
  if (outtime>1000.0){
    outtime = 0.0;
    //cout<<"Time="<<time<<" Sum= "<<sum<<endl;
    f<<"timestep "<<Num<< "  "<<NDsize+2<<endl;
    for(int i=0;i<Num;i++){
      f<<setw(15)<<Z[i]<<setw(15)<<n_neu[i];
      for(int j=0;j<NDsize;j++){
        f<<setw(15)<<Grain[j].man[i];
        }
      f<<endl;
    }
  }
  // updating desorption and adsorption
  for (int i=0;i<NDsize;i++){
  for (int j=1;j<Num-1;j++){
    tneu = Grain[i].inv[j].d1*n_neu[j]+Grain[i].inv[j].d3*Grain[i].man[j];
    tman = Grain[i].inv[j].d2*n_neu[j]+Grain[i].inv[j].d4*Grain[i].man[j];
    //Grain[i].man[j] = tman;
    Grain[i].man[j] += n_neu[j]-tneu;
    n_neu[j] = tneu;
    }
    Grain[i].man[Num-1] = Grain[i].man[Num-2];
    Grain[i].man[0] = Grain[i].man[1];
  }
  // boundary condition
  n_neu[Num-1] = n_neu[Num-2];
  n_neu[0] = n_neu[1];

  if (dif>0.0){ 
  // Convection and diffusion
  for(int j=0;j<Num-1;j++){
   //adflux[j]  = (n_neu[j]*Vz[j]+n_neu[j+1]*Vz[j+1])/2.0; 
   adflux[j]  = n_neu[j+1]*Vz[j+1]/2.0; 
   difflux[j] = -g1[j]*(n_neu[j+1]/N_h[j+1]-n_neu[j]/N_h[j])/dZ;
  }
  adflux[0] = 0.0;
  adflux[Num-2] =0.0;
  difflux[0] = 0.0;
  difflux[Num-2] =0.0;

  for(int j=1;j<Num-1;j++){
    n_neu[j] += dt*(-adflux[j]-difflux[j] +adflux[j-1] +difflux[j-1])/dZ;
  }

  // boundary condition
  n_neu[Num-1] = n_neu[Num-2];
  n_neu[0] = n_neu[1];

  // n_man should follow the same adv/diffu equation as dust
  for (int i=0;i<NDsize;i++){
      for(int j=0;j<Num-1;j++){
         //adflux[j] = (Grain[i].man[j]*Grain[i].Vz[j]+Grain[i].man[j+1]*Grain[i].Vz[j+1])/2.0; 
         adflux[j] = Grain[i].man[j+1]*Grain[i].Vz[j+1]; 
         difflux[j]= -Grain[i].g1[j]*(Grain[i].man[j+1]/N_h[j+1]-Grain[i].man[j]/N_h[j])/dZ;
       }
       adflux[0] = 0.0;
       adflux[Num-2] = 0.0;
       difflux[0] = 0.0;
       difflux[Num-2] =0.0;
      for(int j=1;j<Num-1;j++){
        Grain[i].man[j] += dt*(-adflux[j]-difflux[j] + adflux[j-1] +difflux[j-1])/dZ;
      }
        Grain[i].man[Num-1] = Grain[i].man[Num-2];
        Grain[i].man[0] = Grain[i].man[1];
  }
  }
   
  // enforce positive value and small
  /*
  for (int j=0;j<Num;j++){
     n_neu[j]= max(0.0,n_neu[j]);
     n_neu[j]= min(sum0,n_neu[j]);
   for (int i=0;i<NDsize;i++){
     if(Grain[i].man[j]<0.0){
        n_neu[j] += Grain[i].man[j];
        Grain[i].man[j]= 0.0;
     }else if(Grain[i].man[j]>sum0) {
       n_neu[j] -= sum0-Grain[i].man[j];
       Grain[i].man[j] = sum0;
    }
  }
  }
  */
  // enforce conservation law
  sum = 0.0;
  for (int j=1;j<Num-1;j++){
   sum += n_neu[j];
   for (int i=0;i<NDsize;i++){
      sum += Grain[i].man[j];
    }
  }
  //Grain[NDsize-1].man[1] += (sum0-sum);
  //n_neu[0] += (sum0-sum);
}
/*   */

f<<"final "<<Num<< "  "<<NDsize+2<<endl;
for(int i=0;i<Num;i++){
  f<<setw(15)<<Z[i]<<setw(15)<<n_neu[i];
  for(int j=0;j<NDsize;j++){
    f<<setw(15)<<Grain[j].man[i];
    }
  f<<endl;
}
//cout<<"Evolution succeed!"<<endl;
return 0;
}
