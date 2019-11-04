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
 *         Author:  Francesco Lovascio 
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



real GradDDotGrad(real* P,real* rho, real CS, int i, int j, int k, int pitch, int stride){
	real cs = CS;
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
	FM=1.-(P[llxm]/(CS*CS*rho[llxm]));
	FP=1.-(P[llxp]/(CS*CS*rho[llxp]));
	real x=dx(FP,FM,DX)*dx(P[llxp],P[llxm],DX);
	GdG+=x;
#endif
#ifdef Y
	int llyp=lyp;
	int llym=lym;
	real DYM=ymed(j)-ymed(j-1);
	real DY=ymed(j+1)-ymed(j);
	FM=1.-(P[llym]/(CS*CS*rho[llym]));
	FP=1.-(P[llyp]/(CS*CS*rho[llyp]));
	F=1.-(P[ll]/(CS*CS*rho[ll]));
	real y=dq(FP,F,FM,DY,DYM)*dq(P[llyp],P[ll],P[llym],DY,DYM);
#   ifdef CYLINDRICAL
    y=y/(ymed(j)*ymed(j));
#   endif
	GdG+=y;
#endif
#ifdef Z
	int llzp=lzp;
	int llzm=lzm;
	real DZM=zmed(k)-zmed(k-1);
	real DZ=zmed(k+1)-zmed(k);
	FM=1.-(P[llzm]/(CS*CS*rho[llzm]));
	FP=1.-(P[llzp]/(CS*CS*rho[llzp]));
	F=1.-(P[ll]/(CS*CS*rho[ll]));
	real z=dq(FP,F,FM,DZ,DZM)*dq(P[llzp],P[ll],P[llzm],DZ,DZM);
	GdG+=z;
#endif
	return GdG;
}

real Lap(real* FF, int i, int j, int k, int pitch, int stride){
	//<EXTERNAL>
	real cs = CSCONST;
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
