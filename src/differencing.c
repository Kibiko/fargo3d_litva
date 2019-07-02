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



//#ifdef STABLE
real GradDDotGrad(real* P,real* rho, real CS, int i, int j, int k, int pitch, int stride){
	real cs = CS;
	//real ts = TSCONST;
	////<\EXTERNAL>
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
	real x=dx(FP,FM,DX)*dx(P[llxp],P[llxm],DX)/(H_X(i,j,k)*H_X(i,j,k));
	//(1/(H_X(i,j,k)*H_X(i,j,k)))*((P[llxp]-P[llxm])*(FP-FM)/(4.*DX*DX));
	GdG+=x;
	//	printf("x=%f",x);
#endif
#ifdef Y
	int llyp=lyp;
	int llym=lym;
	real DYM=ymed(j)-ymed(j-1);
	real DY=ymed(j+1)-ymed(j);
	FM=1.-(P[llym]/(CS*CS*rho[llym]));
	FP=1.-(P[llyp]/(CS*CS*rho[llyp]));
	F=1.-(P[ll]/(CS*CS*rho[ll]));
	real y=dq(FP,F,FM,DY,DYM)*dq(P[llyp],P[ll],P[llym],DY,DYM)/(H_Y(i,j,k)*H_Y(i,j,k));
	//(1/(H_Y(i,j,k)*H_Y(i,j,k)))*(DYM/(DY-DYM))*(DY*(P[ll]-P[llym])-DYM*(P[llyp]-P[ll]))*(DYM/(DY-DYM))*(DY*(F-FM)-DYM*(FP-F));
	//	printf("y=%f",y);
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
	real z=dq(FP,F,FM,DZ,DZM)*dq(P[llzp],P[ll],P[llzm],DZ,DZM);//(H_Z(i,j,k)*H_Z(i,j,k));
	//(1/(H_Z(i,j,k)*H_Z(i,j,k)))*(DZM/(DZ-DZM))*(DZ*(P[ll]-P[llzm])-DZM*(P[llzp]-P[ll]))*(DZM/(DZ-DZM))*(DZ*(F-FM)-DZM*(FP-F));
	//	printf("z=%f",z);
	GdG+=z;
#endif
	//printf(" GdG[%d]=%f ",ll,GdG);
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
		//((H_Y(i,j,k)*H_Z(i,j,k)/H_X(i,j,k))*((FF[llxp]+FF[llxm]-2*FF[ll])/(DX*DX)))+
		//(((FF[llxp]-FF[llxm])/(2.*DX)))*(((H_Y(iixp,j,k)*H_Z(iixp,j,k)/H_X(iixp,j,k))-(H_Y(iixm,j,k)*H_Z(iixm,j,k)/H_X(iixm,j,k)))/2.*DX);
	//printf("x=%f",x);
	lap+=x;
#endif
#ifdef Y
	int llyp=lyp;
	int llym=lym;
	real DYM=ymed(j)-ymed(j-1);
	real DY=ymed(j+1)-ymed(j);
	real y=dq2(FF[llyp],FF[ll],FF[llym],DY,DYM);
		//((H_X(i,j,k)*H_Z(i,j,k)/H_Y(i,j,k))*(2./(DY*DY*DYM+DYM*DYM*DY))*((DYM*FF[llyp])+(DY*FF[llym])-((DY+DYM)*FF[ll])))+
		//(DYM/(DY-DYM))*(DY*(FF[ll]-FF[llym])-DYM*(FF[llyp]-FF[ll]))*(DYM/(DY-DYM))*(DY*((H_X(i,j,k)*H_Z(i,j,k)/H_Y(i,j,k))-((H_X(i,j-1,k)*H_Z(i,j-1,k)/H_Y(i,j-1,k))))-DYM*((H_X(i,j+1,k)*H_Z(i,j+1,k)/H_Y(i,j+1,k))-(H_X(i,j,k)*H_Z(i,j,k)/H_Y(i,j,k))));
	//printf("y=%f",y);
	lap+=y;
#endif
#ifdef Z
	int llzp=lzp;
	int llzm=lzm;
	real DZM=zmed(k)-zmed(k-1);
	real DZ=zmed(k+1)-zmed(k);
	real z=dq2(FF[llzp],FF[ll],FF[llzm],DZ,DZM);
		//(((H_X(i,j,k)*H_Y(i,j,k)/H_Z(i,j,k))*2./(DZ*DZ*DZM+DZM*DZM*DZ))*((DZM*FF[llzp])+(DZ*FF[llzm])-((DZ+DZM)*FF[ll])))+
		//(DZM/(DZ-DZM))*(DZ*(FF[ll]-FF[llzm])-DZM*(FF[llzp]-FF[ll]))*(DZM/(DZ-DZM))*(DZ*((H_X(i,j,k)*H_Y(i,j,k)/H_Z(i,j,k))-((H_X(i,j,k-1)*H_Y(i,j,k-1)/H_Z(i,j,k-1))))-DZM*((H_X(i,j,k+1)*H_Y(i,j,k+1)/H_Z(i,j,k+1))-(H_X(i,j,k)*H_Y(i,j,k)/H_Z(i,j,k))));
	//printf("z[%d]=%f ",ll,z);
	lap+=z;
#endif
//	lap*=(1./H_X(i,j,k)*H_Y(i,j,k)*H_Z(i,j,k));
	//printf(" lap[%d]=%f ",ll,lap);
	return lap;
}
#endif
