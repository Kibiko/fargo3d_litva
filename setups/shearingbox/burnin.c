#include "fargo3d.h"


real vortexVZ(real y, real z, real cY, real cZ, real a, real b){
	real v = shearVZinit(y,z);
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
	real v = shearVYinit(y,z);
	real ZZ=z+cZ;
	real YY=y+cY;
	if((ZZ*ZZ/a)+(YY*YY/b)<1.0){
		real q = a/b;
		real w = 1.5*(1.0+q*q)/(q*(q-1));
		v=  w*z/(1+q*q);
	}
	return v;
}

real vortexRho(real y, real z, real cZ, real cY, real a, real b){
   	real v = shearVYinit(y,z);
	real ZZ=z+cZ;
	real YY=y+cY;
	if((ZZ*ZZ/a)+(YY*YY/b)<1.0){
		
	}
	return v;

}

void burnin(real dt) {
  unsigned int size_z=Nz+NGHZ;	
  unsigned int size_y=Ny+NGHY;
  real* vx = Vx->field_cpu;
  real* vy = Vy->field_cpu;
  real* rho= Density->field_cpu;
	
  real maskwidth=KILLZONE;
  real rho0=KILLZONE_RHO_0;
  real e0=KILLZONE_E_0;
  //printf("maskwidth=%f",maskwidth);
  i = j = k = 0;
  for(k=0; k<Nz + 2*NGHZ; k++) {
	for(j=0;j<Ny + 2*NGHY; j++) {
		for(i=0;i<Nx; i++) {
			unsigned int	ll = l;//i+pitch*j+stride*k;
			//printf("%d\n",ll);
			e[ll]		 =e0;
			y1[ll]	 =e0;
			y2[ll]	 =e0;
			y3[ll]	 =e0;
			y4[ll]	 =e0;
			rho[ll]	 =(xmed(i)*xmed(i)/a + ymed(j)*ymed(j)/b)<1.0?4.0*rho0:rho0;
			vy[ll]	 =vortexVYinit(ymed(j),xmed(i),0.0,0.0,a,b);
			vx[ll]	 =vortexVZinit(ymed(j),xmed(i),0.0,0.0,a,b);
		    }
	    }
    }
}



