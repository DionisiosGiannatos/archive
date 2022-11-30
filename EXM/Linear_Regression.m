%One of my earliest university projects utilizing MATLAB/Octave, a linear regression script used in a
%physics experiment
clear;clc
i=0;
X(10)=0;
Y(10)=0;
X=[1.0232 1.0367 1.0396 1.0489 1.0521 1.0547 1.0624 1.0679 1.0963 1.0988]
Y=[3368.32 3362.8 3362.9 3363.7 3363.8 3363.9 3364.8 3364.8 3365.2 3365.6]

Xmean=mean(X);
Ymean=mean(Y);
Xsquared=X.^2;
Xmeansqrd=Xmean^2;
XY=X.*Y;
SumXsqrd=sum(Xsquared);
SumXY=sum(XY);
SSxy=SumXY-(10*Xmean*Ymean);
SSxx=sum(Xsquared)-10*(Xmeansqrd);
%
a=SSxy/SSxx;
b=Ymean-a*Xmean;
