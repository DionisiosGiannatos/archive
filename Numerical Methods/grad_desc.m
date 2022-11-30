#A Gradient Descent algorithm I made for practice, as part of the Numerical Methods class.
#Goal of gradient descent is the location of a minimum
function output = grad_desc (f, xinit,yinit,a=0.3)
hx=sqrt(eps)*abs(xinit);
hy=sqrt(eps)*abs(yinit);
dfdx=@(x,y,hx) (f(x+hx,y)-f(x-hx,y))/(2*hx);
dfdy=@(x,y,hy) (f(x,y+hy)-f(x,y-hy))/(2*hy);
itr=0;
uold=[xinit,yinit];

x(itr+1)=xinit-a*dfdx(xinit,yinit,hx);
y(itr+1)=yinit-a*dfdy(xinit,yinit,hy);

unew=[x(itr+1),y(itr+1)];
chk=norm((unew-uold)./unew,2);
while (chk>10^(-5) && itr<200)
  itr=itr+1;
  uold=[x(itr),y(itr)];
  hx=sqrt(eps)*abs(x(itr));
  hy=sqrt(eps)*abs(y(itr));
  
  x(itr+1)=x(itr)-a*dfdx(x(itr),y(itr),hx);
  y(itr+1)=y(itr)-a*dfdy(x(itr),y(itr),hy);
  
  unew=[x(itr+1),y(itr+1)];
  chk=norm(unew-uold,2);
endwhile
if (itr==200)
  disp("No solution found")
endif
disp("Output: Solution X, Solution Y, Alpha Factor, Number of Iterations")
output=[unew,a,itr]
endfunction
