function Volume = virial(P,T,w,Pc,Tc)
  #Used for the Thermodynamics class (3rd semester), this script accepted 5 arguments: Pressure, Time, Pitzer Acentric Factor
  #w, represetative of the greek omega, critical pressure and critial temperature
  #The original script was lost as octave didn't save the greek characters.
  
  #Reduced Pressure-Temperature Calculations
  Tr = T/Tc;
  Pr = P/Pc;
  R = 8.314; #?? J/K*mol
  
  #Calculation of B for the Virial Expansion.
  
  f0 = 0.1445-(0.33/Tr)-(0.1385/Tr^2)-(0.0121/Tr^3)-(0.000607/Tr^8);
  f1 = 0.0637 + (0.331/Tr^2) - (0.423/Tr^3) - (0.008/Tr^8);
  B = (R*Tc/Pc)*(f0-w*f1); #Eccentricity Factor
  
  #The script below attempts to find the roots of the 2nd degree equation in an attempt to calculate the
  #gas' volume.
  
  #The equation took the form Pv^2 - vRT - BRT, where V was the volume of the gas.
  #It was assumed that Pc<Tc/2, otherwise the truncated Virial equation would give bad results.
  #That being said, the inequality mentioned above is never actually checked in the script.
  
  a=P;
  b= (-1)*R*T;
  c= (-1)*B*R*T;
  coeff = [a,b,c];
  Volume= roots (coeff)
endfunction
