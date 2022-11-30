function output = f30_bonus (K,Mr)
#Calculation of mean speed through the integration of function F(x) 
#Input: Temperature in K, Molecular Weight in g/mol.
#Since this calculates an integral, we substitute infinity with a high number, in this case 2000
#uf(u)=(4/sqrt(pi))*(u/u_mp)^(3)*exp(-(u/u_mp)^(2));
Mr_new=Mr*10^(-3);
u_mp=sqrt(2*8.3141*K/Mr_new);
n=0;
A_new=1;
A_old=2;
while abs((A_new/A_old)-1)>=10^(-8)
 n=n+1;
 A_old=A_new;
 A_tot=0;
 a=0;
 b=4000;
 l=(b-a)/(n);
for i=1:n;
 F_1=erf((a+l*(i-1))/u_mp)-(2/sqrt(pi))*((a+l*(i-1))/u_mp)*exp(-((a+l*(i-1))/u_mp)^2);
 F_2=erf((a+l*i)/u_mp)-(2/sqrt(pi))*((a+l*i)/u_mp)*exp(-((a+l*i)/u_mp)^2);
 A_area=(F_1+F_2)*(l/2);
 A_tot=A_area+A_tot;
endfor
 A_new=A_tot;
 Area_calc=A_new;
endwhile
 u_mean=1*b-Area_calc
 disp("Output: The mean speeed (in m/s) is")
 output=[u_mean]
endfunction
