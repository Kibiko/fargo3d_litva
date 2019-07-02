#LETS CALCULATE STUFF!!
using DelimitedFiles

#FREE PARAMETERS
Rho0=50.0
Rho1=0.01
L=10.0
k=(2.0*pi/L)
ts=1.0
Cs=1.0
f_d=0.5

#DEPENDENT PARAMETERS
F_Omeg()=0.5*((im*(Cs^2)*ts*f_d*k*k)-sqrt(complex(-1*(Cs^4)*(ts^2)*(f_d^2)*(k^4)+4*(k^2)*(Cs^2)*(1-f_d))))
F_P1(Omeg::Complex, P0::Real)=((im*P0)/(im*Omeg+((Cs^2)*ts*k*((1-Rho0)/((Cs^2)*Rho0)))))*(Rho1*Omeg/Rho0)
F_P0()=(Cs^2)*(1-f_d)*(Rho0)
F_V1(Omeg::Complex, P1::Complex)=(k/(Rho0*Omeg))*P1

#PROGRAM LOGIC
P0=F_P0()
Omeg=F_Omeg()
P1=F_P1(Omeg, P0)
V1=F_V1(Omeg,P1)

println(P0)
println(P1)
println(Omeg)
println(V1)

ALL=["Rho0" Rho0;
		 "Rho1" Rho1;
		 "v1" abs(V1);
		 "v1Phase" angle(V1);
		 "P0" P0;
		 "P1" abs(P1);
		 "P1Phase" angle(P1);
		 "k" k;
		 "RE(w)" real(Omeg);
		 "IM(w)" imag(Omeg);
		 "ts" ts;
		 "Cs" Cs;
		 "fd" f_d;
		 "L" L]


#FILE LOGIC
writedlm("WAVE.dat", ALL, "\t")
