//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
#define LEFT  0
#define RIGHT 1
#define DOWN  2
#define UP    3
//<\INCLUDES>

void boundary_zmax_0_cpu () {

//<USER_DEFINED>
INPUT(Density);
INPUT(Energy);
INPUT(Vz);
OUTPUT(Density);
OUTPUT(Energy);
OUTPUT(Vz);
//<\USER_DEFINED>

//<INTERNAL>
  int __attribute__((unused))i;
  int __attribute__((unused))j;
  int __attribute__((unused))k;
  int __attribute__((unused))jact;
  int __attribute__((unused))jgh;
  int __attribute__((unused))kact;
  int __attribute__((unused))kgh;
  int lgh;
  int lghs;
  int lact;
  int lacts;
  int lacts_null;
//<\INTERNAL>

//<EXTERNAL>
  real* rho = Density->field_cpu;
  real* e = Energy->field_cpu;
	real* y1 = Y1->field_cpu;
	real* y2 = Y2->field_cpu;
	real* y3 = Y3->field_cpu;
	real* y4 = Y4->field_cpu;
  real* vz = Vz->field_cpu;
	real* vy = Vy->field_cpu;
  int size_x = Nx+2*NGHX;
  int size_y = Ny+2*NGHY;
  int size_z = NGHZ;
  int nx = Nx;
  int ny = Ny;
  int nz = Nz;
  int nghy = NGHY;
  int nghz = NGHZ;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  real dx = Dx;
//<\EXTERNAL>

//<CONSTANT>
// real xmin(Nx+2*NGHX+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
//<\CONSTANT>

//<MAIN_LOOP>

  i = j = k = 0;

#ifdef Z
  for(k=nz+size_z; k<nz+2*size_z; k++) {
#endif
#ifdef Y
    for(j=0; j<size_y; j++) {
#endif
#ifdef X
      for(i=0; i<size_x; i++) {
#endif
//<#>

	//lgh = i + j*pitch + (nz+nghz+k)*stride;
	//lghs = i + j*pitch + (nz+nghz+1+k)*stride;
	//lact = i + j*pitch + (nz+nghz-1-k)*stride;
	//lacts = i + j*pitch + (nz+nghz-1-k)*stride;
	//lacts_null = i + j*pitch + (nz+nghz)*stride;
	//kgh = (nz+nghz+k);
	//kact = (nz+nghz-1-k);
	//jgh = (ny+nghy+j);

	real rho0=5.0;
	real e0=2.5;
	e[l]		 =e0;
	y1[l]	 =e0;
	y2[l]	 =e0;
	y3[l]	 =e0;
	y4[l]	 =e0;
	rho[l]	 =rho0;
	vy[l]	 =shearVY(ymed(j),zmed(k));
	vz[l]	 =shearVZ(ymed(j),zmed(k));
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
