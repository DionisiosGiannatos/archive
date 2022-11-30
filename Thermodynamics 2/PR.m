#Peng-Robinson EoS Solver
#Calculation of the gas' volume using the Peng Robinson-EoS
function V = PR(P,T,w,Pc,Tc)
  #Calculation of M
  d0 =0.37464;
  d1 =1.54226;
  d2 =(-1)*0.26992;
  R = 8.314;
  m = d0+d1*w+d2*w^2;
  
  #Calculation of Reduced Pressure-Temperature (Pr,Tr)
  Tr = T/Tc;
  Pr = P/Pc;
  aT = (1+m*(1-Tr^0.5))^2;
  
  #Parameters a0,b0
  a0 = 0.45724;
  b0 = 0.07780;
  
  #Calculation of ac,b
  ac = a0*((R*Tc)^2)/Pc;  
  b = b0*(R*Tc)/Pc;
  a=ac*aT;
  
  #Calculation of Peng-Robinson coefficients. the equation solved takes the form of ax^3+bx^2+cx+d=0
  alpha=P;
  beta=(P*b-R*T);
  c=-3*P*(b^2)-2*b*R*T+a;
  delta=P*(b^3)+(b^2)*R*T-a*b;
  
  cubic = [alpha,beta,c,delta];
  
  #Euresh rizwn PR
  
  V=roots(cubic);
endfunction
