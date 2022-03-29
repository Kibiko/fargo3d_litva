#include "fargo3d.h"

void _init_stockholm() {

  INPUT(Density);
  OUTPUT2D(Density0);
#ifdef ADIABATIC
  INPUT(Energy);
  OUTPUT2D(Energy0);
#endif
#ifdef DUSTY
  INPUT(Y1);
  OUTPUT2D(Y1o);
  INPUT(Y2);
  OUTPUT2D(Y2o);
  INPUT(Y3);
  OUTPUT2D(Y3o);
  INPUT(Y4);
  OUTPUT2D(Y4o);
#endif
#ifdef X
  INPUT(Vx);
  OUTPUT2D(Vx0);
#endif
#ifdef Y
  INPUT(Vy);
  OUTPUT2D(Vy0);
#endif
#ifdef Z
  INPUT(Vz);
  OUTPUT2D(Vz0);
#endif

  boolean error_density = TRUE;
  boolean error_vx      = TRUE;
  boolean error_vy      = TRUE;
  boolean error_vz      = TRUE;
  boolean error_energy  = TRUE;

  char outputname[MAXLINELENGTH];

  if ((Restart == YES) || (Restart_Full == YES)) {
    sprintf(outputname,"%s0_2d.dat",Density->name);
    error_density = Read2D(Density0, outputname, OUTPUTDIR, GHOSTINC);
#ifdef X
    sprintf(outputname,"%s0_2d.dat",Vx->name);
    error_vx = Read2D(Vx0, outputname, OUTPUTDIR, GHOSTINC);
#endif
#ifdef Y
    sprintf(outputname,"%s0_2d.dat",Vy->name);
    error_vy = Read2D(Vy0, outputname, OUTPUTDIR, GHOSTINC);
#endif
#ifdef Z
    sprintf(outputname,"%s0_2d.dat",Vz->name);
    error_vz = Read2D(Vz0, outputname, OUTPUTDIR, GHOSTINC);
#endif
#ifdef ADIABATIC
    sprintf(outputname,"%s0_2d.dat",Energy->name);
    error_energy = Read2D(Energy0, outputname, OUTPUTDIR, GHOSTINC);
#endif
  }

  int i,j,k;
  
  i = j = k = 0;
  
#ifdef X
  real* vx  = Vx->field_cpu;
  real* vx0 = Vx0->field_cpu;
#endif
#ifdef Y
  real* vy  = Vy->field_cpu;
  real* vy0 = Vy0->field_cpu;
#endif
#ifdef Z
  real* vz  = Vz->field_cpu;
  real* vz0 = Vz0->field_cpu;
#endif
#ifdef ADIABATIC
  real* e    = Energy->field_cpu;
  real* e0   = Energy0->field_cpu;
#endif
#ifdef DUSTY
  real* y1   = Y1->field_cpu;
  real* y1o  = Y1o->field_cpu;
  real* y2   = Y2->field_cpu;
  real* y2o  = Y2o->field_cpu;
  real* y3   = Y3->field_cpu;
  real* y3o  = Y3o->field_cpu;
  real* y4   = Y4->field_cpu;
  real* y4o  = Y4o->field_cpu;
#endif
  real* rho  = Density->field_cpu;
  real* rho0 = Density0->field_cpu;
  
#ifdef Z
  for (k=0; k<Nz+2*NGHZ; k++) {
#endif
#ifdef Y
    for (j=0; j<Ny+2*NGHY; j++) {
#endif
#ifdef ADIABATIC
      if (error_energy)
	e0[l2D]   = e[l];
#endif
#ifdef DUSTY
      if (error_energy)
	y1o[l2D]  = y1[l];
      	y2o[l2D]  = y2[l];
	y3o[l2D]  = y3[l];
	y4o[l2D]  = y4[l];
#endif
#ifdef X
      if (error_vx)
	vx0[l2D]  = vx[l];
#endif
#ifdef Y
      if (error_vy)
	vy0[l2D]  = vy[l];
#endif
#ifdef Z
      if (error_vz)
	vz0[l2D]  = vz[l];
#endif
      if (error_density)
	rho0[l2D] = rho[l];
#ifdef Y
    }
#endif
#ifdef Z
  }
#endif

  sprintf(outputname,"%s0_2d.dat",Density->name);
  Write2D(Density0, outputname, OUTPUTDIR, GHOSTINC);
#ifdef X
  sprintf(outputname,"%s0_2d.dat",Vx->name);
  Write2D(Vx0, outputname, OUTPUTDIR, GHOSTINC);
#endif
#ifdef Y
  sprintf(outputname,"%s0_2d.dat",Vy->name);
  Write2D(Vy0, outputname, OUTPUTDIR, GHOSTINC);
#endif
#ifdef Z
  sprintf(outputname,"%s0_2d.dat",Vz->name);
  Write2D(Vz0, outputname, OUTPUTDIR, GHOSTINC);
#endif
#ifdef ADIABATIC
  sprintf(outputname,"%s0_2d.dat",Energy->name);
  Write2D(Energy0, outputname, OUTPUTDIR, GHOSTINC);
#endif
#ifdef DUSTY
  sprintf(outputname,"%s0_2d.dat",Y1->name);
  Write2D(Y1o, outputname, OUTPUTDIR, GHOSTINC);
  sprintf(outputname,"%s0_2d.dat",Y2->name);
  Write2D(Y2o, outputname, OUTPUTDIR, GHOSTINC);
  sprintf(outputname,"%s0_2d.dat",Y3->name);
  Write2D(Y3o, outputname, OUTPUTDIR, GHOSTINC);
  sprintf(outputname,"%s0_2d.dat",Y4->name);
  Write2D(Y4o, outputname, OUTPUTDIR, GHOSTINC);
#endif
}

void init_stockholm() {

  static boolean init = TRUE;

  if (init) MULTIFLUID(_init_stockholm());
  
  init = FALSE;

}
