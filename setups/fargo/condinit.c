#include "fargo3d.h"

void Init() {
  
  OUTPUT(Density);
  OUTPUT(Energy);
  OUTPUT(Vx);
  OUTPUT(Vy);

  int i,j,k;
  real r, omega;
  real soundspeed;
  real vmix;
  
  real *vphi = Vx->field_cpu;
  real *vr   = Vy->field_cpu;
  real *rho  = Density->field_cpu;
  real *tsvar = Tsvar->field_cpu;
  
#ifdef ADIABATIC
  real *e   = Energy->field_cpu;
#endif
#ifdef ISOTHERMAL
  real *cs   = Energy->field_cpu;
#endif

  i = j = k = 0;
  
  for (j=0; j<Ny+2*NGHY; j++) {
    for (i=0; i<Nx+2*NGHX; i++) {
      
      r = Ymed(j);
      omega = sqrt(G*MSTAR/r/r/r);
      
      rho[l] = SIGMA0*pow(r/R0,-SIGMASLOPE)*(1.0+NOISE*(drand48()-.5));
      soundspeed  = ASPECTRATIO*pow(r/R0,FLARINGINDEX)*omega*r;

#ifdef ISOTHERMAL
      cs[l] = soundspeed;
#endif
#ifdef ADIABATIC
	#ifndef DUSTY
//	      e[l] = pow(soundspeed,2)*rho[l]/(GAMMA-1.0);
	#endif
#endif
#if defined DUSTY && ADIABATIC
	real *CS = Lics->field_cpu; //18/11
	
//	CS[l] = 0.05;
	CS[l] = soundspeed;
	tsvar[l] = (TSVARIABLE/omega)*(1.-DUSTRATIO);		
	
//	double bump = BUMPTEST;
//	double bumpratio = DUSTRATIO+bump*DUSTRATIO;
//	double cs = CS[l];
//	double rhoold = rho[l];
//	rho[l]+= r>0.8&&r<0.9?bump*DUSTRATIO*rho[l]:0.0;
//	e[l] = r>0.8&&r<0.9?cs*cs*(rho[l]-(1+bump)*DUSTRATIO*rhoold):cs*cs*(1.0-DUSTRATIO)*rhoold;
//	vphi[l] = r>0.8&&r<0.9?sqrt(omega*omega*r*r+cs*cs*(-SIGMASLOPE)*(1.0-bumpratio)/(1.0+bump*DUSTRATIO)):sqrt(omega*omega*r*r+(-SIGMASLOPE)*cs*cs*(1.0-DUSTRATIO)); //CS CONST

	e[l] = rho[l]*CS[l]*CS[l]*(1.0 - DUSTRATIO);
	vphi[l] = omega*r*sqrt(1.0+(1.0-DUSTRATIO)*(2.0*FLARINGINDEX - 1.0 - SIGMASLOPE)*pow(ASPECTRATIO,2.0)*pow(r/R0,2.0*FLARINGINDEX));


	vphi[l] -= OMEGAFRAME*r;
	vr[l] = 0.0;
	
#endif
#ifndef DUSTY      
      vphi[l] = omega*r*sqrt(1.0+pow(ASPECTRATIO,2)*pow(r/R0,2*FLARINGINDEX)*(2.0*FLARINGINDEX - 1.0 - SIGMASLOPE));
      vphi[l] -= OMEGAFRAME*r;
      vphi[l] *= (1.+ASPECTRATIO*NOISE*(drand48()-.5)); //16/12 m

	
      vr[l]    = soundspeed*NOISE*(drand48()-.5);

#endif
    }
  } 
}

void CondInit() {
   Fluids[0] = CreateFluid("gas",GAS);
   SelectFluid(0);
   Init();
}
