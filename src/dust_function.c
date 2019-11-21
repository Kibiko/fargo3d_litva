/*
 * =====================================================================================
 *
 *       Filename:  dust_function.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/20/2018 11:29:19 PM
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



//real C_dust(real yp, real y, real ym, real fdp, real fd, real fdm,real h){
	//<EXTERNAL>
//	real cs = CSCONST; //19/11 m
	//real * CS = LICs->field_cpu; //19/11 a
//	real ts = TSCONST;
	//<\EXTERNAL>

//	real C=0.;
//#ifdef LINEARDIFF
//	C=ts*cs*cs*dxx(yp,y,ym,h)*fd;
//#else
//	C=ts*cs*cs*(dxx(yp,y,ym,h)*fd + dx(fdp,fdm,h)*dx(yp,ym,h));
//#endif
//	C=ts*cs*cs*(dq2(yp,y,ym,h,h)*fd + dq(fdp,fd,fdm,h,h)*dq(yp,y,ym,h,h));
//	return C;
//}

real Cd(real* P_current, real* P, real* rho, real* CS, int i, int j, int k){
	//real cs = CSCONST; //18/11 m2
	real ts = TSCONST;
	int pitch=Pitch_cpu;
	int stride=Stride_cpu;
	real C=0.;
	//C=Lap(P,i,j,k,pitch,stride);
	C=ts*CS[l]*CS[l]*(Lap(P_current,i,j,k,pitch,stride)*(1.-(P_current[l]/(CS[l]*CS[l]*rho[l])))+(GradDDotGrad(P_current,rho,CS,i,j,k,pitch,stride))); //18/11 m2
//	printf("i= %d, lap= %e gdg= %e \n",k,Lap(P_current,i,j,k,pitch,stride)*(1.-(P_current[l]/(cs*cs*rho[l]))), GradDDotGrad(P_current,rho,cs,i,j,k,pitch,stride));
	return C;
}

real gradlncs(real * P_current, real * vy, double flaringindex, int j)
{
	int pitch = Pitch_cpu;
	int stride = Stride_cpu;
	real gradlncs = 0. ;
	gradlncs = -2*(flaringindex-0.5)*P_current[l]*(vy[llym]+vy[llyp])/(2*ymed(j));

	return gradlncs;
}//21/11 a	

#endif
