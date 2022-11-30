function Output = Q_3_Virial (Tc,Pc,Vc,w,y,T,P)
  #This is the calculation of properties of a binary mixture between two gasses, methane and propane
  #Thermodynamic Data is inputed in a (2,1) array
  
  R=8.31446;
  
  for i=1:2;
    #Calculation of the B11, B22, Zc and Z of every material
    Tr(i)=T/Tc(i);
    f0(i)=0.1445-(0.33/Tr(i))-(0.1385/(Tr(i)^2))-(0.0121/(Tr(i)^3))-(0.000607/(Tr(i)^8));
    f1(i)=0.0637+(0.331/Tr(i))-(0.423/Tr(i)^3)-(0.008/(Tr(i)^8));
    B(i)=(R*Tc(i)/Pc(i))*(f0(i)+w(i)*f1(i));
    z(i)=1-(B(i)/(R*T))*P;
    zc(i)=(Pc(i)*Vc(i))/(R*Tc(i));
  endfor
  #Application of combining rules for the two gasses. Methane and propane are considered to be similar
  #thus a l_ij parameter is not required
  Vc_ij=(((Vc(1)^(1/3))+(Vc(2)^(1/3)))/2)^3;
  Tc_ij=(Tc(1)*Tc(2))^(1/2); 
  zc_ij=(zc(1)+zc(2))/2;
  Pc_ij=(zc_ij*R*Tc_ij)/(Vc_ij);
  Tr_ij=T/Tc_ij;
  w_ij=(w(1)+w(2))/2;
  
  #calculating B_12
  f0_ij=0.1445-(0.33/Tr_ij)-(0.1385/(Tr_ij^2))-(0.0121/(Tr_ij^3))-(0.000607/(Tr_ij^8));
  f1_ij=0.0637+(0.331/Tr_ij)-(0.423/Tr_ij^3)-(0.008/(Tr_ij^8));
  B_ij=(R*Tc_ij/Pc_ij)*(f0_ij+w_ij*f1_ij);
  
  #Calculating the observed B
  B_mix=B(1)*(y(1)^2)+(2*y(1)*y(2)*B_ij)+B(2)*(y(2)^2);
  
  #Calculating the Virial parameters (2nd Degree equation)
  alpha=P;
  beta=(-1)*R*T;
  Gamma=(-1)*R*T*B_mix;
  coef=[alpha,beta,Gamma];
  
  #Solution of the Virial equation. because the substance in question is a gas, only the highest root is
  #chosen. Sometimes, roots produced imaginary roots, so only the real part of such roots was chosen
  Vm_mix=roots(coef);
  Vm_actual=real(max(Vm_mix));
  
  #calculating the compressibility factor of the mix through its denfinition
  z_mix=(P*Vm_actual)/(R*T);
  
  #Calculation of the fugacity for materials 1(methane) k 2(propane)
  ln_phi_1=(2/Vm_actual)*(y(1)*B(1)+y(2)*B_ij)-log(z_mix);
  phi_1=exp(ln_phi_1);
  
  ln_phi_2=(2/Vm_actual)*(y(1)*B_ij+y(2)*B(2))-log(z_mix);
  phi_2=exp(ln_phi_2);
  
  
  Output=[Vm_actual,phi_1,phi_2];
  
endfunction
