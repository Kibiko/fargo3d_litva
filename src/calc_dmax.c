/*
 * =====================================================================================
 *
 *       Filename:  calc_dmax.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/20/2018 06:05:14 PM
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

//<INCLUDES>
#include "fargo3d.h"
//#include "./../setups/dust/dust_functions.h"
//</INCLUDES>
#ifdef DUSTY
real calcDmax_cpu(){

	//<USER_DEFINED>
	INPUT(Energy);
	INPUT(Density);
	//INPUT(Mmx);
	//OUTPUT(Mmx);
	//<\USER_DEFINED>

	//<EXTERNAL>
	//double cs = CSCONST; //19/11 m
	real * CS = Lics->field_cpu; //19/11 a
	real * P = Energy->field_cpu;
	real * dens = Density->field_cpu;
	int pitch  = Pitch_cpu;
	int stride = Stride_cpu;
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
	i = j = k = 0;

	//RKL_CASES
	real D_MAX=0.0;
#ifdef Z
	for (k=0; k<size_z; k++) {
#endif
#ifdef Y
		for (j=0; j<size_y; j++) {
#endif
#ifdef X
			for (i=0; i<size_x; i++ ) {
#endif
				//<#>
				ll = l;
				real f_d= 1-(P[ll]/(CS[ll]*CS[ll]*dens[ll])); //19/11 a
				if(f_d>D_MAX){
					D_MAX=f_d;
				}
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
	return D_MAX;
}

#endif //DUSTY
