#include "fargo3d.h"
#include <stdio.h>
#define FLOAT_PI 3.14159265358979 
#define FN "./setups/dustyWave/WAVE.dat"


void Init() {
  
	int i,j,k;
  real* rho = Density->field_cpu;
  real* e = Energy->field_cpu;
	real* y1 = Y1->field_cpu;
	real* y2 = Y2->field_cpu;
	real* y3 = Y3->field_cpu;
	real* y4 = Y4->field_cpu;

#ifdef DWAVE
   double Data[14];/*  =
	 {50.0,
	  0.01,
	  7.520123344669839e-5,
	  -2.115852999896501,
	  25.0,
	  0.002658765106187651,
		0.8017266382316853,
	  0.6283185307179586,
	  1.,
	  0.,
		1.0,
		1.0,
		0.5,
		10.0
		};*/

 	FILE* fp = fopen(FN,"r");

  double v;
  int n;
  for(n=0;n<14;n++){
	fscanf(fp,"%*s");
	fscanf(fp,"%lf",&v);
	//printf("%lf\n",v);
	Data[n]=v;
  }
  fclose(fp);
	double cs = Data[11];
	double _fd = Data[12];
  int index;
#ifdef Z
  real* v1 = Vz->field_cpu;
  int dim1 = Nz + 2*NGHZ;
#define Q1 (Zmed(i) - ZMIN)/(ZMAX - ZMIN)
#define E1 (Zmin(i) - ZMIN)/(ZMAX - ZMIN)
#define N1 Nz
#endif
#ifdef Y
  real* v1 = Vy->field_cpu;
  int dim1 = Ny + 2*NGHY;
#define Q1 (Ymed(i) - YMIN)/(YMAX - YMIN)
#endif
#ifdef X
  real* v1 = Vx->field_cpu;
  int dim1 = Nx;
#define Q1 (Xmed(i) - XMIN)/(XMAX - XMIN)
#define E1 (Xmin(i) - XMIN)/(XMAX - XMIN)
#endif
    for (i = 0; i<dim1; i++) {
		//printf("%f\n",Data[13]*Q1);
    e[i]   =Data[4]-(Data[5]*cos((Q1)*Data[7]*Data[13] +Data[6]));
    rho[i] =Data[0]-(Data[1]*cos(Q1*Data[7]*Data[13]));
   	v1[i]  =Data[2]*cos((E1)*Data[7]*Data[13]+Data[3]);
	}
#endif //DWAVE

#ifdef DSHOCK
   double Rho0=10;
   double CS2=1;
	 double fd=0.5;
   double P0=-CS2*(fd-1)*Rho0;
#ifdef Z
  real* v1 = Vz->field_cpu; 
  int dim1 = Nz + 2*NGHZ;
#define Q1 (Zmed(i) - ZMIN)/(ZMAX - ZMIN)
#endif
#ifdef Y
  real* v1 = Vy->field_cpu;
  int dim1 = Ny + 2*NGHY;
#define Q1 (Ymed(i) - YMIN)/(YMAX - YMIN)
#endif
#ifdef X
  real* v1 = Vx->field_cpu;
  int dim1 = Nx;
#define Q1 (Xmed(i) - XMIN)/(XMAX - XMIN)
#endif
  
  for (i = 0; i<dim1; i++) {
    
		e[i] = (Q1<0.5?0.5*P0:P0);
		rho[i] = (Q1<0.5?0.4*Rho0:Rho0);

		
	//	e[i]   = P0;
    //rho[i] = Rho0;
    v1[i]  = 0.0;
  //  if (Q1 > 0.4) {
	//		rho[i] = 0.25*Rho0;
	//		e[i] = 0.5*P0;
	//	}
//		y1[i]=0;
//		y2[i]=0;
//		y3[i]=0;
//		y4[i]=0;
	}
#endif //DSHOCK
#ifdef RWAVE
   double Rho0=10;
   double CS2=2;
   double P0=CS2*0.9*Rho0;
#ifdef Z
  real* v1 = Vz->field_cpu; 
  int dim1 = Nz + 2*NGHZ;
#define Q1 (Zmed(i) - ZMIN)/(ZMAX - ZMIN)
#endif
#ifdef Y
  real* v1 = Vy->field_cpu;
  int dim1 = Ny + 2*NGHY;
#define Q1 (Ymed(i) - YMIN)/(YMAX - YMIN)
#endif
#ifdef X
  real* v1 = Vx->field_cpu;
  int dim1 = Nx;
#define Q1 (Xmed(i) - XMIN)/(XMAX - XMIN)
#endif
  
  for (i = 0; i<dim1; i++) {
    e[i]   = P0;
    rho[i] = Rho0;
    v1[i]  = 0.0;
    if (Q1 > 0.5) {
		v1[i]=1.0;
	}
	}
#endif //RWAVE
#ifdef SHOCK
#ifdef Z
  real* v1 = Vz->field_cpu;
  int dim1 = Nz + 2*NGHZ;
#define Q1 (Zmed(i) - ZMIN)/(ZMAX - ZMIN)
#endif
#ifdef Y
  real* v1 = Vy->field_cpu;
  int dim1 = Ny + 2*NGHY;
#define Q1 (Ymed(i) - YMIN)/(YMAX - YMIN)
#endif
#ifdef X
  real* v1 = Vx->field_cpu;
  int dim1 = Nx;
#define Q1 (Xmed(i) - XMIN)/(XMAX - XMIN)
#endif
  
  for (i = 0; i<dim1; i++) {
    e[i]   = 1.0;
    rho[i] = 1.0;
    v1[i]  = 0.0;
    if (Q1 > 0.5) {
      rho[i] = 0.125;
	  e[i] = 0.125;
    }
  }
#endif //SHOCK
}

void CondInit() {
   Fluids[0] = CreateFluid("gas",GAS);
   SelectFluid(0);
   Init();
}
