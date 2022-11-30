function phi = PR_fugacity(P,vPR,T,w,Pc,Tc)
  #Calculation of Fugacity using the Peng-Robinson Equation of State (EoS)
  #Calculation of M
  d0 =0.37464;
  d1 =1.54226;
  d2 =(-1)*0.26992;
  R = 8.314;
  m = d0+d1*w+d2*w^2;
  
  #Calculation of reduced Pressure-Temperature
  Tr = T/Tc;
  Pr = P/Pc;
  aT = (1+m*(1-Tr^0.5))^2;
  
  #a0,b0 Parameters
  a0 = 0.45724;
  b0 = 0.07780;
  
  #Calculation of ac,b Parameters
  ac = a0*((R*Tc)^2)/Pc;  
  b = b0*(R*Tc)/Pc;
  a=ac*aT;
  
  #Calculation of A,B,zV parameters
  A=(a*P)/((R*T)^2);
  B=(b*P)/(R*T);
  zgas=(P*vPR)/(R*T);
  
  #Calculation of Fugacity
  Lnphi=(zgas-1)-(log(zgas-B))-((A)/(B*2*sqrt(2)))*log((zgas+(1+sqrt(2))*B)/(zgas+(1-sqrt(2))*B));
  phi=exp(Lnphi);
endfunction
