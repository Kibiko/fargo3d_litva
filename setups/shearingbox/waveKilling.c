#include "fargo3d.h"

real shearVY(real y, real x){
	real v = 0.0;
	return v;
}

real shearVX(real y, real x){
	real v = -1.5*y;
	return v;
}

real maskfunction(real distance_ratio){
	real bet=0.008;
	return bet*sin(0.5*3.141592653597932*distance_ratio)*sin(0.5*3.141592653597932*distance_ratio);
}

void waveKiller(real dt) {
//	printf("entering wavekiller\n"); 
	int i,j,k;
  real* rho = Density->field_cpu;
  real* e = Energy->field_cpu;
	unsigned int size_z=Nz+NGHZ;	
	unsigned int size_y=Ny+NGHY;
  real* vx = Vx->field_cpu;
  real* vy = Vy->field_cpu;
	
	real maskwidth=KILLZONE;
	real rho0=KILLZONE_RHO_0;
	real e0=KILLZONE_E_0;
	//printf("maskwidth=%f",maskwidth);
	i = j = k = 0;
	for(k=NGHZ; k<size_z; k++) {
		for(j=NGHY;j<size_y; j++) {
			for(i=0;i<Nx; i++) {

				unsigned int	ll = l;
				if( (Ymed(j)-YMIN-(Ymed(j)-Ymin(j)))<maskwidth){
					real Q=(Ymed(j)-YMIN+(Ymed(j)-Ymin(j)))/maskwidth;
					vy[ll] -= (maskfunction(Q))*(vy[ll]-shearVY(Ymed(j),Xmed(i)));	//maskfunction(Q)*vy[ll]+
									//(1.0-maskfunction(Q))*shearVY(Ymin(j),Xmin(i));
					vx[ll] -= (maskfunction(Q))*(vx[ll]-shearVX(Ymed(j),Xmed(i)));	//maskfunction(Q)*vx[ll]+
									//(1.0-maskfunction(Q))*shearVX(Ymin(j),Xmin(i));
					e[ll] -= 	(maskfunction(Q))*(e[ll]-e0);//maskfunction(Q)*e[ll]+
									//(1.0-maskfunction(Q))*e0;
					rho[ll] -= 	(maskfunction(Q))*(rho[ll]-rho0);//maskfunction(Q)*rho[ll]+
									//(1.0-maskfunction(Q))*rho0;
				}
				else if(Ymed(j)>(YMAX-maskwidth+(Ymed(j)-Ymin(j)))){
					real Q=(YMAX-Ymed(j)+(Ymed(j)-Ymin(j)))/maskwidth;
					vy[ll] -= 	(maskfunction(Q))*(vy[ll]-shearVY(Ymed(j),Xmed(i)));//maskfunction(Q)*vy[ll]+
									//(1.0-	maskfunction(Q))*shearVY(Ymin(j),Xmin(i));
					vx[ll] -= 	(maskfunction(Q))*(vx[ll]-shearVX(Ymed(j),Xmed(i)));//maskfunction(Q)*vx[ll]+
									//(1.0-maskfunction(Q))*shearVX(Ymin(j),Xmin(i));
					e[ll] -= 	(maskfunction(Q))*(e[ll]-e0);//maskfunction(Q)*e[ll]+
									//(1.0-maskfunction(Q))*e0;
					rho[ll] -= (maskfunction(Q))*(rho[ll]-rho0);//maskfunction(Q)*rho[ll]+
									//(1.0-maskfunction(Q))*rho0;
				}
			}
		}
	}
}



