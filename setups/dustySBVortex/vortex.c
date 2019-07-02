#include "fargo3d.h"
#define FLOAT_PI 3.14159265358979 

real shearVY(real y, real z){
	real v = 0.0;
	return v;
}

real shearVZ(real y, real z){
	real v = -1.5*y;
	return v;
}

real vortexVZ(real y, real z){
	real v = shearVZ(y,z);
	return v;
}

real vortexVY(real y, real z){
	real v = shearVY(y,z);
	return v;
}

