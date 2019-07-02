#Dusty vortex testing script

mpirun -np=6 ./fargo3d setups/shearingbox/shear_TS05.par
mpirun -np=6 ./fargo3d setups/shearingbox/shear_TS01.par
mpirun -np=6 ./fargo3d setups/shearingbox/shear_TS001.par
