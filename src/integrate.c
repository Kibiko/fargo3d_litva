/*
 * =====================================================================================
 *
 *       Filename:  integrate.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/15/2018 16:04:01
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>
#ifdef DUSTY
//<INCLUDES>
#include "fargo3d.h"
//#include "./../setups/dust/dust_functions.h"
//</INCLUDES>

void integrate_cpu(real dt, real * F1, real * F2, real * Fnew, int J, int S){

	//<USER_DEFINED>
	//<\USER_DEFINED>

	//<EXTERNAL>
	//double safety = SAFETYFACTOR;
	int size_x = XIP;
	int size_y = Ny+2 * NGHY;
	int size_z = Nz+2 * NGHZ;
	real RhoDust = RHODUST;
	real SurfDust = SURFDUST;
	real TS_CONST = TSCONST;
	real * rho = Density->field_cpu;
	real * CS = Lics->field_cpu; //18/11 a
	real * F0 = Energy->field_cpu;
	real * store = Y4->field_cpu;
	real * vy = Vy_temp->field_cpu; //21/11 a
	//<\EXTERNAL>

	//<INTERNAL>
	int i;
	int j;
	int k;
	int ll;
	double divP;
	//<\INTERNAL>

	//<CONSTANT>
	// real xmin(Ny+2*NGHX+1);
	// real ymin(Ny+2*NGHY+1);
	// real zmin(Nz+2*NGHZ+1);
	// real GAMMA(1);
	//<\CONSTANT>


	//<MAIN_LOOP>
//printf("dt=%e \n",dt);
#ifdef Z
	for (k= NGHZ-1; k<size_z-1; k++) {
#endif
#ifdef Y
		for (j=NGHY-1; j<size_y-1; j++) {
#endif
#ifdef X
			for (i=0; i<size_x; i++ ) {
#endif
				
				ll = l;	
				Fnew[ll]=0.;
#	ifdef RKL1
				Fnew[ll] = ((2.*J-1)/J)*F1[ll]+ 
				((1.-J)/J)*F2[ll]+ 
			  //((1.0-muj(J)-nuj(J))*F0[ll])+ 
				((4.*J-2.)/(J*S*S+J*S))*dt*Cd(F1,F1,rho,CS,i,j,k); //18/11 a
				//store[ll] = Fnew[ll];
#	else
				Fnew[ll] = (muj(J)*F1[ll])+(nuj(J)*F2[ll])+((1.0-muj(J)-nuj(J))*F0[ll])+(mu1j(J,S)*dt*Cd(F1,F1,rho,CS,i,j,k))+(nu1j(J,S)*dt*Cd(F1,F0,rho,CS,i,j,k));
#	endif			
#ifdef X
			}
#endif
#ifdef Y
		}
#endif
#ifdef Z
	}

#endif
}

void integrate1_cpu(real dt, real * F0, real * Fnew, int S){

	//<USER_DEFINED>
	//<\USER_DEFINED>

	//<EXTERNAL>
	//double safety = SAFETYFACTOR;
	int size_x = XIP;
	int size_y = Ny+2 * NGHY;
	int size_z = Nz+2 * NGHZ;
	real RhoDust = RHODUST;
	real SurfDust = SURFDUST;
	real TS_CONST = TSCONST;
	real * rho = Density->field_cpu;
	real * CS = Lics->field_cpu; //18/11 a
	//<\EXTERNAL>

	//<INTERNAL>
	int i;
	int j;
	int k;
	int ll;
	double divP;
	real *y0;
	real *y1;
	real *ym;
	real *ymm;
	//<\INTERNAL>

	//<CONSTANT>
	// real xmin(Ny+2*NGHX+1);
	// real ymin(Ny+2*NGHY+1);
	// real zmin(Nz+2*NGHZ+1);
	// real GAMMA(1);
	//<\CONSTANT>


	//<MAIN_LOOP>
#ifdef Z
	for (k=NGHZ-1 ; k<size_z-NGHZ; k++) {
#endif
#ifdef Y
		for (j=NGHY-1; j<size_y-NGHY; j++) {
#endif
#ifdef X
			for (i=0; i<size_x; i++ ) {
#endif
				ll = l;
				Fnew[ll]=0.;
				Fnew[ll]+= F0[ll]+mu1j(1,S)*dt*Cd(F0,F0,rho,CS,i,j,k);
#ifdef X
			}
#endif
#ifdef Y
		}
#endif
#ifdef Z
	}
#endif
}

void RK2_cpu(real dt){
	//<USER_DEFINED>
	//<\USER_DEFINED>

	//<EXTERNAL>
	//double safety = SAFETYFACTOR;
	int size_x = Nx;
	int size_y = Ny+2 * NGHY;
	int size_z = Nz+2 * NGHZ;
	real RhoDust = RHODUST;
	real SurfDust = SURFDUST;
	real TS_CONST = TSCONST;
  real * rho = Density->field_cpu;
	real * F0 = Energy->field_cpu;
	real * F1 = Y1->field_cpu;
	real * CS = Lics->field_cpu; //18/11 a
	//<\EXTERNAL>

	//<INTERNAL>
	int i;
	int j;
	int k;
	int ll;
	double divP;
	int comunicate;  
	//<\INTERNAL>

	//<CONSTANT>
	// real xmin(Ny+2*NGHX+1);
	// real ymin(Ny+2*NGHY+1);
	// real zmin(Nz+2*NGHZ+1);
	// real GAMMA(1);
	//<\CONSTANT>


	//<MAIN_LOOP>
//printf("dt=%e \n",dt);
#ifdef Z
	for (k=NGHZ-1; k<size_z-NGHZ; k++) {
#endif
#ifdef Y
		for (j=NGHY-1; j<size_y-NGHY; j++) {
#endif
#ifdef X
			for (i=0; i<size_x; i++ ) {
#endif
				//<#>
				ll = l;
				F1[ll]=F0[ll]+0.5*dt*Cd(F0,F0,rho,CS,i,j,k); //18/11 a
				//<\#>
#ifdef X
			}
#endif
#ifdef Y
		}
#endif
	current_simulation_time+=0.5*dt;
	comunicate = Y1_COMM;
	FARGO_SAFE(FillGhosts(comunicate));

#ifdef Z
	}
#endif
#ifdef Z
	for (k=NGHZ-1; k<size_z-NGHZ; k++) {
#endif
#ifdef Y
		for (j=NGHY-1; j<size_y-NGHY; j++) {
#endif
#ifdef X
			for (i=0; i<size_x; i++ ) {
#endif
				//<#>
				ll = l;
				F0[ll]=F0[ll]+dt*Cd(F1,F1,rho,CS,i,j,k); //18/11
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
	current_simulation_time+=0.5*dt;
	comunicate = ENERGY;
	FARGO_SAFE(FillGhosts(comunicate));

}

void Assign_cpu1(real * E, real * Enew){

	//<USER_DEFINED>
	//<\USER_DEFINED>

	//<EXTERNAL>
	//double safety = SAFETYFACTOR;
	int size_x = Nx;
	int size_y = Ny+2 * NGHY;
	int size_z = Nz+2 * NGHZ;
	real RhoDust = RHODUST;
	real SurfDust = SURFDUST;
	real TS_CONST = TSCONST;


	//<\EXTERNAL>

	//<INTERNAL>
	int i;
	int j;
	int k;
	int ll;
	double divP;
	real *y0;
	real *y1;
	real *ym;
	real *ymm;
	//<\INTERNAL>

	//<CONSTANT>
	// real xmin(Ny+2*NGHX+1);
	// real ymin(Ny+2*NGHY+1);
	// real zmin(Nz+2*NGHZ+1);
	// real GAMMA(1);
	//<\CONSTANT>


	//<MAIN_LOOP>
#ifdef Z
	for (k=NGHZ; k<size_z-NGHZ; k++) {
#endif
#ifdef Y
		for (j=NGHY; j<size_y-NGHY; j++) {
#endif
#ifdef X
			for (i=0; i<size_x; i++ ) {
#endif
				//<#>
				ll = l;
				//printf("DeltaE[%d]=%e \n",ll,Enew[ll]-E[ll]); 
				E[ll] = Enew[ll]; //11/12 m
				
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
	//printf("\n");
}

#endif
