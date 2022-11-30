function [output] = Series_2(T,Tc,Pc,w,Radius,Grams)
  #This was a script used to solve the 2nd series of exercises for the class of Thermodynamics 2
  #and one of my first university scripts to utilize an iterative method.
  #Judging from the script, this was used for vapor-liquid equilibrium calculations
  
  #For a random pressure
  P=rand();
  #Calculate M
  d0 =0.37464;
  d1 =1.54226;
  d2 =(-1)*0.26992;
  R = 8.314;
  m = d0+(d1*w)+d2*(w^2);
  Vmin=10000000; Vmax=(-1);
  f_liquid=2;
  f_vapor=1;
  
  #Calculation of total Moles and Volume
  Mr=58.124 #g*(mol^(-1))
  Vol_Cont=(4/3)*pi*(Radius^3);
  mol_total=(Grams)/(Mr);
  
  #The iterative method. The goal is to have the calculated fugacity of liquid and vapor phases
  #be as close as possible to each other.
  while ((f_liquid/f_vapor)-1)>=10^(-6) 
   #Calculation of Reduced Pressure and Reduced Temperature
   Pnew=P*(f_liquid/f_vapor);
   P=Pnew;
   Tr = T/Tc;
   Pr = P/Pc;
   aT = (1+m*(1-Tr^0.5))^2;
  
   #a0,b0 parameters
   a0 = 0.45724;
   b0 = 0.07780;
  
   #Calculation of ac,b
   ac = a0*((R*Tc)^2)/Pc;  
   b = b0*(R*Tc)/Pc;
   a=ac*aT;
   
   #Calculation of A,B
   A=(a*P)/((R*T)^2);
   B=(b*P)/(R*T);
   
   #Solution of the Peng-Robinson equation. The equation has the form of ax^3+bx^2+cx+d=0
   alpha=P;
   beta=(P*b-R*T);
   c=-3*P*(b^2)-2*b*R*T+a;
   delta=P*(b^3)+(b^2)*R*T-a*b;
  
   cubic = [alpha,beta,c,delta];
  
   #The actual solution of PR EoS, used to calculate the volume
   V=roots(cubic);
   V2=real(V);
   
   #A convoluted script, whose purpose is to keep the min and max roots while discarding the middle solution
   #Somehow, the min and max commands eluded me.
   for i=1:3
     Vtemp=V2(i,1);
     if (Vtemp>Vmax) && (Vtemp>0);
       Vmax=Vtemp;
     endif
     if (Vtemp<Vmin) && (Vtemp>0);
       Vmin=Vtemp;
     endif
   endfor
   #keeping the min and max solutions.
   VmL=Vmin;
   VmV=Vmax;
   Vmax=0;
   Vmin=10000000; #A "reset" of the values of Vmax, Vmin.
   zliq=(P*VmL)/(R*T);
   zvap=(P*VmV)/(R*T);
   
   #calculation of fugacity coefficients
   phi_vap=exp((zvap-1)-log(zvap-B)-(A/(2*B*sqrt(2)))*log((zvap+(1+sqrt(2))*B)/(zvap+(1-sqrt(2))*B)));
   phi_liq=exp((zliq-1)-log(zliq-B)-(A/(2*B*sqrt(2)))*log((zliq+(1+sqrt(2))*B)/(zliq+(1-sqrt(2))*B)));
   
   #calculation of fugacities
   f_liquid=P*phi_liq;
   f_vapor=P*phi_vap;
  endwhile
  Ps=P;
  #calculation of moles in each phase
  mol_liq=(Vol_Cont-(mol_total*VmV))/(VmL-VmV);
  mol_g=mol_total-mol_liq;
  
  #calculation of mass in each phase.
  kg_g=mol_g*Mr*10^(-3)
  kg_liq=mol_liq*Mr*10^(-3)
  
  output=[Ps,VmV,VmL,kg_g,kg_liq];
endfunction