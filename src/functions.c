//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//#include "./../setups/dust/dust_functions.h"
//</INCLUDES>

#ifdef DUSTY

real bj(int J){
	if(J<3){
		return 0.33333333333333;
	}
	else{
		return (1.0*J*J+J-2)/(2.0*J*J-2.0*J);
	}
}

real ws(int s){
	return 4.0/(s*s+s-2);
}

real aj(int J){
	return 1.0 - bj(J);
}

real muj(int J){
	return ((2.0*J-1.0)/(1.0*J))*(bj(J)/bj(J-1));
}

real nuj(int J){
	return -1.0*((1.0*J-1.0)/(1.0*J))*(bj(J)/bj(J-2));
}

real mu1j(int J, int s){
	return muj(J)*ws(s);
}

real nu1j(int J, int s){
	return -1.0*(1.0-bj(J))*mu1j(J,s);
}

#	ifdef STABLE
void make_D_cpu(){

	//<USER_DEFINED>
	INPUT(Energy);
	INPUT(Density);
	//INPUT(Mmx);
	//OUTPUT(Mmx);
	//<\USER_DEFINED>

	//<EXTERNAL>
	real * D = Y4->field_cpu;
	real * p = Energy->field_cpu;
	real * rho=Density->field_cpu;
	int pitch  = Pitch_cpu;
	int stride = Stride_cpu;
	int size_x = Nx;
	int size_y = Ny+2 * NGHY;
	int size_z = Nz+2 * NGHZ;

	//real cs = CSCONST; //19/11 m
	real * CS = LICs->field_cpu; //19/11 a
	real ts = TSCONST;
	//<\EXTERNAL>

	//<INTERNAL>
	int i;
	int j;
	int k;
	int ll;
	//<\INTERNAL>

	//<CONSTANT>
	// real xmin(Ny+2*NGHX+1);
	// real ymin(Ny+2*NGHY+1);
	// real zmin(Nz+2*NGHZ+1);
	// real GAMMA(1);
	//<\CONSTANT>


	//<MAIN_LOOP>

#		ifdef Z
	for (k=0; k<size_z; k++) {
#		endif
#		ifdef Y
		for (j=0; j<size_y; j++) {
#		endif
#		ifdef X
			for (i=0; i<size_x; i++ ) {
#		endif
				//<#>
				ll = l;
				D[ll]=1-(p[ll]/(CS[ll]*CS[ll]*rho[ll])); //19/11 m
				//<\#>
#		ifdef X
			}
#		endif
#		ifdef Y
		}
#		endif
#		ifdef Z
	}
#		endif
}

#	endif
#endif

