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
	real *CS = LICs->field_cpu; //18/11
	CS[l] = soundspeed;
	e[l] = rho[l]*pow(CS[l],2.)*(1.-DUSTRATIO); //18/11 stores Pressure in Energy field
	vmix = omega*r*sqrt(1.0+(1.0-DUSTRATIO)*(2.0*FLARINGINDEX - 1.0 - SIGMASLOPE)*pow(ASPECTRATIO,2.0)*pow(r/R0,2.0*FLARINGINDEX));
	vphi[l] = vmix;
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
