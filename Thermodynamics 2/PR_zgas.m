function zgas = PR_zgas(P,T,w,Pc,Tc)
  #Calculation of M
  d0 =0.37464;
  d1 =1.54226;
  d2 =(-1)*0.26992;
  R = 8.314;
  m = d0+d1*w+d2*w^2;
  
  #Calculation of Tr kai Pr
  Tr = T/Tc;
  Pr = P/Pc;
  aT = (1+m*(1-Tr^0.5))^2;
  
  #a0,b0 Parameters
  a0 = 0.45724;
  b0 = 0.07780;
  
  #Calculation of ac,b
  ac = a0*((R*Tc)^2)/Pc;  
  b = b0*(R*Tc)/Pc;
  a=ac*aT;
  
  #Calculation of A,B,zV
  A=a*P/((R*T)^2);
  B=b*P/(R*T);
  
  #Solving the Peng Robinson EoS to identify the compressibility factor Z
  alpha=1;
  beta=(1-B);
  gama=(A-2*B-3*(B^2));
  delta=(A*B-B^2-B^3);
  
  COEF=[alpha,beta,gama,delta];
  zgas=roots(COEF);
endfunction
