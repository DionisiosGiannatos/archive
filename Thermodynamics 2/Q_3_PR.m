function Output = Q_3_PR(Tc,Pc,Vc,w,y,T,P)
  #The same process as Q_3_Virial, used to solve question 3. The main difference is the usage of the 
  #Peng-Robinson EoS instead of Virial
  
  #Definition of constants
  R=8.31446;
  
  d0=0.37464;
  d1=1.54226;
  d2= -0.26992;
  
  a0=0.45724;
  b0=0.07780;
  
  #Calculation of a,b for every compound. A and B will be later used for the fugacity coefficients.

  for i=1:2
    Tr(i)=T/Tc(i);
    m(i)=d0+d1*w(i)+d2*(w(i)^2);
   
    a_T(i)=(1+m(i)*(1-(Tr(i)^(0.5))))^(2);
    a_c(i)=a0*((R*Tc(i))^2)/Pc(i);
    
    b(i)=b0*(R*Tc(i))/Pc(i);
    a(i)=a_c(i)*a_T(i);
    
    A(i)=a(i)*P/((R*T)^2);
    B(i)=b(i)*P/(R*T);
  endfor
  #Application of combining rules for the calculation of a,b,A,B of the compounds
  
  a_ij=(a(1)*a(2))^(1/2);
  a_mix=a(1)*(y(1)^2)+(2*y(1)*y(2)*a_ij)+a(2)*(y(2)^2);
  
  b_ij=(b(1)+b(2))/2;
  b_mix=b(i)*(y(1)^2)+(2*y(1)*y(2)*b_ij)+b(2)*(y(2)^2);
  
  A_ij=a_ij*P/((R*T)^2);
  B_ij=b_ij*P/(R*T);
  
  alpha=P;
  beta=(P*b_mix-R*T);
  c=-3*P*(b_mix^2)-2*b_mix*R*T+a_mix;
  delta=P*(b_mix^3)+(b_mix^2)*R*T-a_mix*b_mix;
  
  cubic = [alpha,beta,c,delta];
  Vm_solution=roots(cubic);
  Vm_mix=max(real(Vm_solution));
  
  z_mix= (P*Vm_mix)/(R*T);
  
  A_mix=a_mix*P/((R*T)^2);
  B_mix=b_mix*P/(R*T);
  
  ln_phi_1=(B(1)/B_mix)*(z_mix-1)-log(z_mix-B_mix)-(A_mix/(2*sqrt(2)*B_mix))*((2*(y(1)*A(1)+(y(2)*A_ij))/(A_mix))-(B(1)/B_mix))*log((z_mix+(1+sqrt(2))*B_mix)/(z_mix+(1-sqrt(2))*B_mix));
  phi_1=exp(ln_phi_1);
  
  ln_phi_2=(B(2)/B_mix)*(z_mix-1)-log(z_mix-B_mix)-(A_mix/(2*sqrt(2)*B_mix))*((2*(y(1)*A_ij+(y(2)*A(2)))/(A_mix))-(B(2)/B_mix))*log((z_mix+(1+sqrt(2))*B_mix)/(z_mix+(1-sqrt(2))*B_mix));
  phi_2=exp(ln_phi_2);
  
  Output=[Vm_mix, phi_1,phi_2];
endfunction
