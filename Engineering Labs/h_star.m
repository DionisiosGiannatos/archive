function y = h_star(x)
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
#x(1)=C_H*, x(2)=CO_2
a=zeros(2,1);
a(1)=Cad*((Ka1/x(1))+2*(Ka1*Ka2)/(x(1)^2)+3*(Ka1*Ka2*Ka3)/(x(1)^(3)))/(1+(Ka1/(x(1)))+(Ka1*Ka2)/(x(2)^2)+(Ka1*Ka2*Ka3)/(x(1)^3))-x(1);
a(2)=Kb1*Kb2*(x(1)/2)*(sqrt(1+(4*Cbd)/(Kb1*(1+Kb2*x(1))))-1);
y=zeros(2,1);
y(1)=a(1)-a(2);
y(2)=a(1)-x(2);
endfunction
