function y = citric(x)
#Constants
pH=2.5;
Ka1=10^(-3.13);
Ka2=10^(-4.76);
Ka3=10^(-6.4);
H=10^(-pH);
#Equations
#x(1)=Cc, x(2)=Cc1, x(3)=Cc2, x(4)=Cc3, x(5)=CA_0
y=zeros(5,1);
y(1)=x(2)-Ka1*x(1)/H;
y(2)=x(3)-Ka2*x(2)/H;
y(3)=x(4)-Ka3*x(3)/H;
y(4)=x(2)+2*x(3)+3*x(4)-H;
y(5)=x(5)-x(4)-x(3)-x(2)-x(1);
endfunction
