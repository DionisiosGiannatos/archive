function y = tot_eq(x)
#Constants
Ca=1.6488e-02; #Original Citric Concentration
Cb=0.1; #Original Soda Concentration
Qa=5.3; #Citric Acid Flow Rate
Qb=10.3; #Sodium Bicarbonate Flow Rate
Cad=Ca*(Qa/(Qa+Qb));#Ca diluted from mixing with soda
Cbd=Cb*(Qb/(Qa+Qb)); #Cb diluted from mixing with citric acid
Ka1=10^(-3.13);
Ka2=10^(-4.76);
Ka3=10^(-6.4);
Kb1=0.69;
Kb2=10^(6.351);
#Equations
#x(1)=C_H*, x(2)=Cc, x(3)=Cc1, x(4)=Cc2, x(5)=Cc3, x(6)=NaHCO3, x(7)=HCO3-, x(8)=CO2, x(9)=Na+ (concentrations)
y=zeros(9,1);
y(1)=Ka1-x(3)*x(1)/x(2);
y(2)=Ka2-x(4)*x(1)/x(3);
y(3)=Ka3-x(5)*x(1)/x(4);
y(4)=Kb1-x(7)*x(9)/x(6);
y(5)=Kb2-x(8)/(x(7)*x(1));
y(6)=x(2)+x(3)+x(4)+x(5)-Cad;
y(7)=x(6)+x(9)-Cbd;
y(8)=x(6)+x(7)+x(8)-Cbd;
y(9)=x(1)+x(9)-x(3)-2*x(4)-3*x(5)-x(7);
endfunction
