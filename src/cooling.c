
#include "fargo3d.h"

#ifdef DUSTY

void cooling_cpu(real dt){
#ifndef TESTNOCOOLING
  SelectFluid(0);
	//double cs = CSCONST; //19/11 m
	real * CS = LICs->field_cpu; //19/11 a
	real * gradlc = glcs->field_cpu; //10/12 a
	real * y0 = Energy->field_cpu;
	real * y1 = Y1->field_cpu;
	real * y2 = Y2->field_cpu;
	real * y3 = Y3->field_cpu;
	real * y4 = Y4->field_cpu;
	real * vel = Vz->field_cpu;
	real * dens = Density->field_cpu;
	int pitch  = Pitch_cpu;
	int stride = Stride_cpu;
	int size_x = Nx;
	int size_y = Ny+2 * NGHY;
	int size_z = Nz+2 * NGHZ;
	real RhoDust = RHODUST;
	real SurfDust = SURFDUST;
	real TS_CONST = TSCONST;
	int i;
	int j;
	int k;
	int ll;
	double divP;
	
	i = j = k = 0;
	real t_parabX=dt;
	real t_parabY=dt;
	real t_parabZ=dt;

	real DMAX=TS_CONST*CS[ll]*CS[ll]; //19/11 m

#ifdef Z
	real t_parabz=0.1*(zmed(2)-zmed(1))*(zmed(2)-zmed(1))/DMAX;
#ifdef FLOAT
	MPI_Allreduce(&t_parabz, &t_parabZ, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
#else
	MPI_Allreduce(&t_parabz, &t_parabZ, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif
#endif

#ifdef Y
	real t_paraby=0.1*(ymed(2)-ymed(1))*(ymed(2)-ymed(1))/DMAX;
#ifdef FLOAT
	MPI_Allreduce(&t_paraby, &t_parabY, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
#else
	MPI_Allreduce(&t_paraby, &t_parabY, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif
#endif

#ifdef X 
	real t_parabx=0.1*(xmed(2)-xmed(1))*(xmed(2)-xmed(1))/DMAX;
#ifdef FLOAT
	MPI_Allreduce(&t_parabx, &t_parabX, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
#else
	MPI_Allreduce(&t_parabx, &t_parabX, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif
#endif
	
	real t_parab=t_parabX>t_parabY?t_parabY:t_parabX;
	t_parab=t_parab>t_parabZ?t_parabZ:t_parab;

	

	int intersteps=ceil(fabs(0.5*(sqrt(9+16*dt/t_parab)-1)));
	//printf("[intersteps=%d]\n",intersteps);
	//printf("t_parab= %f \n",t_parab);
	//printf("DMAX= %f \n",DMAX);
	double time_initial=current_simulation_time;
	int J=2;
	FARGO_SAFE(Assign_cpu1(y1,y0))
	FARGO_SAFE(Assign_cpu1(y2,y0))
	FARGO_SAFE(Assign_cpu1(y3,y0))
	FARGO_SAFE(Assign_cpu1(y4,y0))
	FARGO_SAFE(integrate1_cpu(dt, y0, y1, intersteps));
	int comunicate=Y1_COMM;
	//comunicate |= ENERGY;
	FARGO_SAFE(FillGhosts(comunicate));
	if(J>1){
		current_simulation_time=time_initial +dt*(4.0/((double)(intersteps*intersteps+intersteps-2)));
		FARGO_SAFE(integrate_cpu(dt, y1, y0, y2, J, intersteps));  
		comunicate=Y2_COMM;
	}
	//comunicate |= ENERGY;
	FARGO_SAFE(FillGhosts(comunicate));
	if(intersteps>2){
		for(J=3;J<intersteps+1;J++){
			current_simulation_time=time_initial +dt*((double)(J*J+J-2)/((double)(intersteps*intersteps+intersteps-2)));
			//printf("J=%d\n",J);
			if((J%3)==0){
				FARGO_SAFE(integrate_cpu(dt, y2, y1, y3, J, intersteps));
				int comunicate=Y3_COMM;
				//		comunicate |= ENERGY;
				FARGO_SAFE(FillGhosts(comunicate));
				//Assign_cpu(y4,y3);
			}
			if((J%3)==1){
				FARGO_SAFE(integrate_cpu(dt, y3, y2, y1, J, intersteps));
				int comunicate=Y1_COMM;
				//		comunicate |= ENERGY;
				FARGO_SAFE(FillGhosts(comunicate));
				//Assign_cpu(y4,y1);
			}
			if((J%3)==2){
				FARGO_SAFE(integrate_cpu(dt, y1, y3, y2, J, intersteps));
				int comunicate=Y2_COMM;
				//		comunicate |= ENERGY;
				FARGO_SAFE(FillGhosts(comunicate));	
				//Assign_cpu(y4,y2);
			}
			//printf("current time %d = %f\n",J,current_simulation_time);
		}
			//printf("\n \n");
		}
	if((J%3)==0){
		Assign_cpu1(y0,y3);
		//Assign_cpu(y0,y3,gradlc);
	}
	if((J%3)==1){
		Assign_cpu1(y0,y1);
		//Assign_cpu(y0,y1,gradlc);
	}
	if((J%3)==2){
		Assign_cpu1(y0,y2);
		//Assign_cpu(y0,y2,gradlc);
	}
	//printf("\n \n");
	current_simulation_time=time_initial+dt;
	comunicate = ENERGY;
	FARGO_SAFE(FillGhosts(comunicate));
#endif
}


void RK2_cooling_cpu(real dt){

	//double cs = CSCONST; //19/11 m
	real * CS = LICs->field_cpu;
	real * y0 = Energy->field_cpu;
	real * dens = Density->field_cpu;
	int pitch  = Pitch_cpu;
	int stride = Stride_cpu;
	int size_x = Nx;
	int size_y = Ny+2 * NGHY;
	int size_z = Nz+2 * NGHZ;
	real RhoDust = RHODUST;
	real SurfDust = SURFDUST;
	real TS_CONST = TSCONST;
	int i;
	int j;
	int k;
	int ll;
	double divP;
	
	i = j = k = 0;
	real DeltaT=0.0;
	real DMAX=TS_CONST*CS[ll]*CS[ll]*calcDmax_cpu(); //19/11 m
	real t_parabX=dt;
	real t_parabY=dt;
	real t_parabZ=dt;
#ifdef Z
	real t_parab1=0.3*(zmed(2)-zmed(1))*(zmed(2)-zmed(1))/DMAX;
#ifdef FLOAT
	MPI_Allreduce(&t_parab1, &t_parabZ, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
#else
	MPI_Allreduce(&t_parab1, &t_parabZ, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif
#endif

#ifdef Y
	real t_parab2=0.3*(ymed(2)-ymed(1))*(ymed(2)-ymed(1))/DMAX;
#ifdef FLOAT
	MPI_Allreduce(&t_parab2, &t_parabY, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
#else
	MPI_Allreduce(&t_parab2, &t_parabY, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif
#endif

#ifdef X 
	real t_parab3=0.3*(xmed(2)-xmed(1))*(xmed(2)-xmed(1))/DMAX;
#ifdef FLOAT
	MPI_Allreduce(&t_parab3, &t_parabX, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
#else
	MPI_Allreduce(&t_parab3, &t_parabX, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif
#endif
	real t_parab=t_parabX>t_parabY?t_parabY:t_parabX;
	t_parab=t_parab>t_parabZ?t_parabZ:t_parab;
	DeltaT+=t_parab;
	
	FARGO_SAFE(RK2_cpu(t_parab));

#ifdef PRINT_STEPS
	printf("dt/tau=%f\ndelta t=%f",t_parab/dt,DeltaT);
#endif
	
	while(DeltaT<dt){

		real DMAX=TS_CONST*CS[ll]*CS[ll]*calcDmax_cpu(); //19/11 m
		t_parabX=dt-DeltaT;
		t_parabY=dt-DeltaT;
		t_parabZ=dt-DeltaT;
#ifdef Z
		real t_parab1=0.3*(zmed(2)-zmed(1))*(zmed(2)-zmed(1))/DMAX;
#ifdef FLOAT
		MPI_Allreduce(&t_parab1, &t_parabZ, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
#else
		MPI_Allreduce(&t_parab1, &t_parabZ, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif
#endif

#ifdef Y
		real t_parab2=0.3*(ymed(2)-ymed(1))*(ymed(2)-ymed(1))/DMAX;
#ifdef FLOAT
		MPI_Allreduce(&t_parab2, &t_parabY, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
#else
		MPI_Allreduce(&t_parab2, &t_parabY, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif
#endif

#ifdef X 
		real t_parab3=0.3*(xmed(2)-xmed(1))*(xmed(2)-xmed(1))/DMAX;
#ifdef FLOAT
		MPI_Allreduce(&t_parab3, &t_parabX, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
#else
		MPI_Allreduce(&t_parab3, &t_parabX, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif
#endif
		real t_parab=t_parabX>t_parabY?t_parabY:t_parabX;
		t_parab=t_parab>t_parabZ?t_parabZ:t_parab;
		DeltaT+=t_parab;
#ifdef PRINT_STEPS
		printf("dt/tau=%f",t_parab/dt);
#endif
		if(DeltaT>=dt){break;}
		FARGO_SAFE(RK2_cpu(t_parab));
	}
	t_parab=fabs(dt-DeltaT);
	FARGO_SAFE(RK2_cpu(t_parab));


}
#endif//DUSTY
