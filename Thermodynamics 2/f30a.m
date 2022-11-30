function output = f30a(x)
#calculation of the error function through approximative methods (riemann sum) from 0 to 1
#x: value that we want to calculate
#n: Natural number: the amount of repetitions that we want.
n=0;
A_new=1;
A_old=2;
while abs((A_new/A_old)-1)>=10^(-7)
 n=n+1;
 A_old=A_new;
 A_tot=0;
 a=0;
 f_a=(2/sqrt((pi)))*e^(0);
 b=x;
 f_b=(2/sqrt(pi))*e^(-b^2);
 l=(b-a)/(n);
for i=1:n;
 f_1=(2/(sqrt(pi)))*e^(-(a+(i-1)*l)^2);
 f_2=(2/sqrt(pi))*e^(-(a+i*l)^2);
 A_area=(f_1+f_2)*(l/2);
 A_tot=A_area+A_tot;
endfor
 A_new=A_tot;
endwhile
Area_calc=A_new;
disp("Output type: Number of Itterations, Area Calculated");
output=[n,Area_calc];
endfunction
