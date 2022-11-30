clear; clc;
%Solution of the Wind Generation heat problem through the usage of Analytical Methods. All errors are
%Expressed in Percentages.
%Constant Initialization
fprintf("Solution of the Wind Turbine problem. \n\nTemperature found in class: \n")
T_real=405.13
D=3; L=6; h=35; e=0.83; s=5.67*10^(-8);
A=(pi*D*L+(2/4)*pi*D^2);
Tsur=20+273.15;
Tinf=25+273.15;
Q=329654.7821;
i=0;
error=1;
%Function Initialization. These functions are used throughout the script.
f=@(T) A*(e*s*(T^(4)-Tsur^(4))+h*(T-Tinf))-Q;
hs=@(T) e*s*(T^(2)+Tsur^(2))*(T+Tsur);
%Calculation of T through Fzero
 fprintf("\nCalculation of Temperature through the fzero function\n")
T_zero=fzero(f,[-500,500])
Error_zero=100*abs(T_real-T_zero)/T_real
%Calculation of T by ignoring Heat Transfer through Radiation
fprintf("\nCalculation of Temperature by ignoring Heat Transfer through Radiation\n")
T_no_rad=Tinf+Q/(A*h)
T(1)=T_no_rad;
Error_no_rad=100*abs(T_real-T_no_rad)/T_real
%Calculation through Repetition. A linearized method is utilized to express the non-linear
  fprintf("\nCalculation of Temperature through a Repetitive Method\n")
  fprintf("Starting Temp: %d\n\n", T(1))
while (error>=10^(-7)) && (i<200);
  i=i+1;
  fprintf("Repetition: %d\n",i)
  const=hs(T(i));
  fprintf("Hs factor: %d\n",const)
  T(i+1)=((Q/A)+const*Tsur+h*Tinf)/(const+h);
  fprintf("Current Temperature: %f\n\n",(T(i+1))) 
  error=abs(T(i+1)-T(i))/T(i+1);
endwhile
  fprintf("Final Temperature and Error\n")
  T_rep=T(i)
  Error_rep=100*abs(T_rep-T_real)/T_real
