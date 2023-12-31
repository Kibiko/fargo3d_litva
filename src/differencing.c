/*
 * =====================================================================================
 *
 *       Filename:  differencing.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/20/2018 06:05:14 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Francesco Lovascio and Kevin Chan 
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

real dx(real ip, real im, real h){
	return (ip-im)/(2*h);
}

real dxx(real ip, real i, real im, real h){
	return (ip+im-2*i)/(h*h);
}

real dq(real fp, real f, real fm, real hi, real hm){
	return (hm*hm*fp-hi*hi*fm+(hi*hi-hm*hm)*f)/(hi*hm*(hi+hm));
}

real dq2(real fp, real f, real fm, real hi, real hm){
	return 2*(hm*fp+hi*fm-f*(hi+hm))/(hi*hm*(hi+hm));
}

real GradDDotGrad(real* P,real* rho, real* CS, real* TS, int i, int j, int k, int pitch, int stride)//18/11 m
	{
	//real cs = CSCONST; //18/11 m
	//real ts = TSCONST;
	real GdG=0.;
	int ll =l;
	real FM;
	real F;
	real FP;
#ifdef X
	int llxp=lxp;
	int llxm=lxm;
	real DX=zone_size_x(j,k);
	FM=TS[llxm]*(1.-(P[llxm]/(CS[llxm]*CS[llxm]*rho[llxm]))); //18/11 m
	FP=TS[llxp]*(1.-(P[llxp]/(CS[llxp]*CS[llxp]*rho[llxp]))); //18/11 m
	real x=dx(FP,FM,DX)*dx(P[llxp],P[llxm],DX);
#ifdef CYLINDRICAL
	x = x/(ymed(j)*ymed(j));
#endif
	GdG+=x;
#endif
#ifdef Y
	int llyp=lyp;
	int llym=lym;
	real DYM=ymed(j)-ymed(j-1);
	real DY=ymed(j+1)-ymed(j);
	FM=TS[llym]*(1.-(P[llym]/(CS[llym]*CS[llym]*rho[llym]))); //18/11 m
	FP=TS[llyp]*(1.-(P[llyp]/(CS[llyp]*CS[llyp]*rho[llyp]))); //18/11 m
	F=TS[ll]*(1.-(P[ll]/(CS[ll]*CS[ll]*rho[ll]))); //18/11 m
	real y=dq(FP,F,FM,DY,DYM)*dq(P[llyp],P[ll],P[llym],DY,DYM);

	GdG+=y;
#endif
#ifdef Z
	int llzp=lzp;
	int llzm=lzm;
	real DZM=zmed(k)-zmed(k-1);
	real DZ=zmed(k+1)-zmed(k);
	FM=TS[llzm]*(1.-(P[llzm]/(CS[llzm]*CS[llzm]*rho[llzm]))); //18/11 m
	FP=TS[llzp]*(1.-(P[llzp]/(CS[llzp]*CS[llzp]*rho[llzp]))); //18/11 m
	F=TS[ll]*(1.-(P[ll]/(CS[ll]*CS[ll]*rho[ll]))); //18/11 m
	real z=dq(FP,F,FM,DZ,DZM)*dq(P[llzp],P[ll],P[llzm],DZ,DZM);
	GdG+=z;
#endif
	return GdG;
}

real Lap(real* FF, int i, int j, int k, int pitch, int stride){
	//<EXTERNAL>
	//real cs = CSCONST; //18/11 m
	real ts = TSCONST;
	//<\EXTERNAL>
	real lap=0.;
	int ll =l;
#ifdef X
	int llxp=lxp;
	int llxm=lxm;
	int iixp=ixp;
	int iixm=ixm;
	real DX=zone_size_x(j,k);
	real x=dxx(FF[llxp],FF[ll],FF[llxm],DX);
#   ifdef CYLINDRICAL
    x=x/(ymed(j)*ymed(j));
#   endif
	lap+=x;
#endif
#ifdef Y
	int llyp=lyp;
	int llym=lym;
	real DYM=ymed(j)-ymed(j-1);
	real DY=ymed(j+1)-ymed(j);
	real y=dq2(FF[llyp],FF[ll],FF[llym],DY,DYM);
#ifdef CYLINDRICAL
	y=y+dq(FF[llyp],FF[ll],FF[llym],DY,DYM)/ymed(j);
#endif
	lap+=y;
#endif
#ifdef Z
	int llzp=lzp;
	int llzm=lzm;
	real DZM=zmed(k)-zmed(k-1);
	real DZ=zmed(k+1)-zmed(k);
	real z=dq2(FF[llzp],FF[ll],FF[llzm],DZ,DZM);
	lap+=z;
#endif
	return lap;
}

#endif
