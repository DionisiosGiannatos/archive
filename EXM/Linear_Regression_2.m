%Another Linear Regression script for another physics experiment. It also calculates the coefficient of
%determination (R^2) of the fit.
%Looking back four years later, it might not be the most efficient of scripts, but it got the work done
clear;clc
i=0;
X(10)=0;
Y(10)=0;
X=[1.0232 1.0367 1.0396 1.0489 1.0521 1.0547 1.0624 1.0679 1.0963 1.0988]
Y=[3368.32 3362.8 3362.9 3363.7 3363.8 3363.9 3364.8 3364.8 3365.2 3365.6]
%Linear Regression
Xsqrd=X.^2;
XY=X.*Y;
Xmean=mean(X);
Ymean=mean(Y);
SumXsqrd=sum(Xsqrd);
SumXY=sum(XY);
SSxyprecursor=(10*Xmean*Ymean);
SSxy=SumXY-(SSxyprecursor); 
SSxx=SumXsqrd-10*(Xmean)^2;
%Calculation of parameters
a=SSxy/SSxx;
disp('a is: '), disp(a)
b=Ymean-a*Xmean;
disp('b is: '), disp(b)
disp(['The form of the equation is: y=', num2str(a),"x+", num2str(b)]) %Displaying the results string
for i=1:10
  Yest(i)=a*X(i)+b;
endfor
%Calculation of the S
for i=1:10
  SSEprecursor(i)=(Y(i)- Yest(i))^2;
endfor
SSE=sum(SSEprecursor);
%Calculation of SST
for i=1:10
  SSTprecursor(i)=(Y(i)-Ymean)^2;
endfor
SST=sum(SSTprecursor);
Rsqrd=1-(SSE/SST) %Calculation of R^2
disp('The values estimated by the liniar regression are: '), disp(Yest)
Ydif=Y-Yest; %
disp('The difference between the values given and the values given by the linear regression are: '), disp(Ydif)