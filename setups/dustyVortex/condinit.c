#include "fargo3d.h"
#define FLOAT_PI 3.14159265358979 

real shearVYinit(real y, real z){
	real v = 0.0;
	return v;
}

real shearVZinit(real y, real z){
	real v = -1.5*y;
	return v;
}

real vortexVZinit(real y, real z, real cY, real cZ, real a, real b){
	real v = shearVZinit(y,z);
	real ZZ=z+cZ;
	real YY=y+cY;
	if((ZZ*ZZ/a)+(YY*YY/b)<1.0){
		real q = a/b;
		real w = 1.5*(1.0+q*q)/(q*(q-1));
	//	v= - q*q*w*y/(1+q*q);
	}
	return v;
}

real vortexVYinit(real y, real z, real cZ, real cY, real a, real b){
	real v = shearVYinit(y,z);
	real ZZ=z+cZ;
	real YY=y+cY;
	if((ZZ*ZZ/a)+(YY*YY/b)<1.0){
		real q = a/b;
		real w = 1.5*(1.0+q*q)/(q*(q-1));
	//	v=  w*z/(1+q*q);
	}
	return v;
}

void Init() {
  
	int i,j,k;
  real* rho = Density->field_cpu;
  real* e = Energy->field_cpu;
	real* y1 = Y1->field_cpu;
	real* y2 = Y2->field_cpu;
	real* y3 = Y3->field_cpu;
	real* y4 = Y4->field_cpu;
	current_simulation_time=0.0;

	real a = 0.05;
	real b = 0.0125;

  real* vz = Vz->field_cpu;
  int dim1 = Nz + 2*NGHZ;
  real* vy = Vy->field_cpu;
	real rho0=5.0;
	real e0=25.0;
	i = j = k = 0;
	for(k=0; k<Nz + 2*NGHZ; k++) {
		for(j=0;j<Ny + 2*NGHY; j++) {
		unsigned int	ll = l;
				e[ll]		 =e0;
				y1[ll]	 =e0;
				y2[ll]	 =e0;
				y3[ll]	 =e0;
				y4[ll]	 =e0;
				rho[ll]	 =rho0;
				vy[l]	 =vortexVYinit(ymed(j),zmed(k),0.0,0.0,a,b);
				vz[l]	 =vortexVZinit(ymed(j),zmed(k),0.0,0.0,a,b);
		}
//		printf("\n");
	}
}

void CondInit() {
   Fluids[0] = CreateFluid("gas",GAS);
   SelectFluid(0);
   Init();
}
