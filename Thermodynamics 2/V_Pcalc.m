function [output]=V_Pcalc(T,Tc,Pc,w)
  #Calculation of Volume for a random Pressure
  #I am not sure how this script is different from Series_2, only that it actually presents the data
  #to the person running the function.
  
  P=rand();
  #Calculation of M
  d0 =0.37464;
  d1 =1.54226;
  d2 =(-1)*0.26992;
  R = 8.314;
  m = d0+(d1*w)+d2*(w^2);
  Vmin=10000000; Vmax=(-1);
  f_liquid=2;
  f_vapor=1;
  k=0;
  #Initialization of the iterative loop
  while ((f_liquid/f_vapor)-1)>=10^(-6) 
   k=k+1;
   #Calculation of Tr kai Pr
   Pnew=P*(f_liquid/f_vapor);
   P=Pnew;
   Tr = T/Tc;
   Pr = P/Pc;
   aT = (1+m*(1-Tr^0.5))^2;
  
   #a0,b0 parameters
   a0 = 0.45724;
   b0 = 0.07780;
  
   #ypologismos ac,b
   ac = a0*((R*Tc)^2)/Pc;  
   b = b0*(R*Tc)/Pc;
   a=ac*aT;
   
   #Calculation of A,B
   A=(a*P)/((R*T)^2);
   B=(b*P)/(R*T);
   
   #Coefficients of the Peng-Robinson equation. The equation is in the form of ax^3+bx^2+cx+d=0
   alpha=P;
   beta=(P*b-R*T);
   c=-3*P*(b^2)-2*b*R*T+a;
   delta=P*(b^3)+(b^2)*R*T-a*b;
  
   cubic = [alpha,beta,c,delta];
  
   #Finding the roots of the PR EoS
   V=roots(cubic);
   V2=real(V);
   
   #Selecting the min kai max root.
   for i=1:3
     Vtemp=V2(i,1);
     if (Vtemp>Vmax) && (Vtemp>0);
       Vmax=Vtemp;
     endif
     if (Vtemp<Vmin) && (Vtemp>0);
       Vmin=Vtemp;
     endif
   endfor
   #VmL and VmV are the volumes of the liquid and vapor phase.
   VmL=Vmin;
   VmV=Vmax;
   Vmax=0;
   Vmin=10000000; #Resetting Min and Max
   zliq=(P*VmL)/(R*T);
   zvap=(P*VmV)/(R*T);
   
   #calculation of fugacity coefficients
   phi_vap=exp((zvap-1)-log(zvap-B)-(A/(2*B*sqrt(2)))*log((zvap+(1+sqrt(2))*B)/(zvap+(1-sqrt(2))*B)));
   phi_liq=exp((zliq-1)-log(zliq-B)-(A/(2*B*sqrt(2)))*log((zliq+(1+sqrt(2))*B)/(zliq+(1-sqrt(2))*B)));
   
   #calculation of the fugacities.
   f_liquid=P*phi_liq;
   f_vapor=P*phi_vap;
  endwhile
  Ps=P;
  output=[VmV,VmL,Ps];
  disp("Number of itterations:"), disp(k);
endfunction
