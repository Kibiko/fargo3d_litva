//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void SubStep3_cpu (real dt) {

//<USER_DEFINED>
  INPUT(Energy);
#ifdef X
  INPUT(Vx_temp);
#endif
#ifdef Y
  INPUT(Vy_temp);
#endif
#ifdef Z
  INPUT(Vz_temp);
#endif
  OUTPUT(Energy);
//<\USER_DEFINED>

//<EXTERNAL>
  real* e   = Energy->field_cpu;
  real * CS = Lics->field_cpu;
#ifdef X
  real* vx  = Vx_temp->field_cpu;
#endif
#ifdef Y
  real* vy  = Vy_temp->field_cpu;
  real* vy_ini = Vy->field_cpu;
#endif
#ifdef Z
  real* vz  = Vz_temp->field_cpu;
#endif
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = XIP; 
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
//<\EXTERNAL>

//<INTERNAL>
  int i; //Variables reserved
  int j; //for the topology
  int k; //of the kernels
  int ll;
#ifdef X
  int llxp;
  int llxm;
#endif
#ifdef Y
  int llyp;
  int llym;
#endif
#ifdef Z
  int llzp;
#endif
  real term;
  real div_v;
  real gradlc; //6/12 a
//<\INTERNAL>
  
//<CONSTANT>
// real GAMMA(1);
// real Sxj(Ny+2*NGHY);
// real Syj(Ny+2*NGHY);
// real Szj(Ny+2*NGHY);
// real Sxk(Nz+2*NGHZ);
// real Syk(Nz+2*NGHZ);
// real Szk(Nz+2*NGHZ);
// real InvVj(Ny+2*NGHY);
//<\CONSTANT>

//<MAIN_LOOP>
  
  i = j = k = 0;
  
#ifdef Z
  for(k=0; k<size_z; k++) {
#endif
#ifdef Y
    for(j=NGHY-1; j<size_y-NGHY; j++) {
#endif
#ifdef X
      for(i=0; i<size_x; i++) {
#endif
//<#>

	ll = l;
#ifdef X
	llxp = lxp;
	llxm = lxm;
#endif
#ifdef Y
	llyp = lyp;
	llym = lym;
#endif
#ifdef Z
	llzp = lzp;
#endif

#ifndef DUSTY
	div_v = 0.0;
#ifdef X
	div_v += (vx[llxp]-vx[ll])*SurfX(j,k);
#endif
#ifdef Y
	div_v += (vy[llyp]*SurfY(j+1,k)-vy[ll]*SurfY(j,k));
#endif
#ifdef Z
	div_v += (vz[llzp]*SurfZ(j,k+1)-vz[ll]*SurfZ(j,k));
#endif
	term = 0.5 * dt * (GAMMA- 1.) * div_v * InvVol(j,k);
	e[ll] *= (1.0-term)/(1.0+term);
#endif

#ifdef DUSTY
#ifndef TESTNOGRAD
//	div_v = 0.0;
//#ifdef X
//        div_v += (vx[llxp]-vx[ll])*SurfX(j,k);
//#endif
//#ifdef Y
//        div_v += (vy[llyp]*SurfY(j+1,k)-vy[ll]*SurfY(j,k));
//#endif
//#ifdef Z
//        div_v += (vz[llzp]*SurfZ(j,k+1)-vz[ll]*SurfZ(j,k));
//#endif
//	term = 0.5 * dt * div_v * InvVol(j,k);
//	e[ll] *= (1.0-term)/(1.0+term);

	gradlc = 0.0;
#ifdef Y
        real r = Ymed(j);
//	gradlc = 1.0/(1.0-(dt*vy[ll]*(- 1.0)/r)); //original 
//	gradlc = 1.0-dt*0.5*(vy[ll]+vy[llyp])/r; //test1 explicit
	gradlc = exp(-0.5*(vy[ll]+vy[llyp])*dt/r); //test2 Chen and Lin 
#endif
	e[ll] *= gradlc; 

#endif
#endif

//<\#>
#ifdef X
      }
#endif
#ifdef Y
    }
#endif
#ifdef Z
  }
#endif
//<\MAIN_LOOP>
}
