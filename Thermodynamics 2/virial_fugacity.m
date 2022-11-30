function phi = virial_fugacity(P,vVirial,T,w,Pc,Tc)
#Calculation of fugacity through the Virial EoS
  #Basic parameters
  Tr = T/Tc;
  Pr = P/Pc;
  R = 8.314;
  
  #Calculation of B through the Tsonopoulos method
  f0 = 0.1445-(0.33/Tr)-(0.1385/Tr^2)-(0.0121/Tr^3)-(0.000607/Tr^8);
  f1 = 0.0637 + (0.331/Tr^2) - (0.423/Tr^3) - (0.008/Tr^8);
  B = (R*Tc/Pc)*(f0-w*f1);
  
  #Calculation of fugacity
  z=1+(B/vVirial);
  Lnphi=(2*B/vVirial)-log(z);
  phi=exp(Lnphi);
endfunction
