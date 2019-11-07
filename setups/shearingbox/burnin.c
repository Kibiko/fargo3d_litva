#include "fargo3d.h"

real vortexVZ(real y, real z, real cY, real cZ, real a, real b){
  real v = shearVX(y,z);
  real ZZ=z+cZ;
  real YY=y+cY;
  if((ZZ*ZZ/a)+(YY*YY/b)<1.0){
    real q = a/b;
    real w = 1.5*(1.0+q*q)/(q*(q-1));
    v= - q*q*w*y/(1+q*q);
  }
  return v;
}

real vortexVY(real y, real z, real cZ, real cY, real a, real b){
  real v = shearVY(y,z);
  real ZZ=z+cZ;
  real YY=y+cY;
  if((ZZ*ZZ/a)+(YY*YY/b)<1.0){
    real q = a/b;
    real w = 1.5*(1.0+q*q)/(q*(q-1));
    v=  w*z/(1+q*q);
  }
  return v;
}

real burninSmoothingFunction(int nt){
  real sf=cos(0.5*M_PI*((double)nt/((double)BURNINSTEPS)));
  sf*=sf;
  return sf;
}

void burnin(int nt) {
  if(nt<BURNINSTEPS){
    unsigned int size_z=Nz+NGHZ;	
    unsigned int size_y=Ny+NGHY;
    real* vx = Vx->field_cpu;
    real* vy = Vy->field_cpu;
    real* rho= Density->field_cpu;
    real maskwidth=KILLZONE;
    real rho0=KILLZONE_RHO_0;
    real e0=KILLZONE_E_0;
    real a=vortexA;
    real b=vortexB;
    int i;
    int j;
    int k;
    double bsf=burninSmoothingFunction(nt); //An attempt to make relaxation profile smooth in time
    i = j = k = 0;
    for(k=0; k<Nz + 2*NGHZ; k++) {
      for(j=0;j<Ny + 2*NGHY; j++) {
        for(i=0;i<Nx; i++) {
          unsigned int	ll = l;//i+pitch*j+stride*k;
          rho[ll]	 =bsf*rho[ll]+(1.0-bsf)*((xmed(i)*xmed(i)/a + ymed(j)*ymed(j)/b)<1.0?4.0*rho0:rho0);
          vy[ll]	 =bsf*vy[ll]+(1.0-bsf)*vortexVY(ymed(j),xmed(i),0.0,0.0,a,b);
          vx[ll]	 =bsf*vx[ll]+(1.0-bsf)*vortexVZ(ymed(j),xmed(i),0.0,0.0,a,b);
        }
      }
    }
  }
}



