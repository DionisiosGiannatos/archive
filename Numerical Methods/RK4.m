#Runge-Kutta 4th degree calculator.
#follows the classic Runge-Kutta method of y(j+1)=y(j)+DT((1/6)k1+(2/6)k2+(2/6)k3+(1/6)k4)
#the problem has the following form: dy/dt=f(t,y) or similar, ex: dy/dx=f(x,y)
#fhandler is the input function, must be of the form f(x,y) or f(t,x)
#only problem is the insertion of y0 into 0, if one wants to calculate the solution at negative t or x.
# Theoretically it would require a lot of checks and a "backwards" R-K
function output=RK4(fhandler,y0,xmin=0,xmax=1,points=100)
  x=linspace(xmin,xmax,points);
  h=x(3)-x(2);
  y(1)=y0;
  for i=1:(points-1)
    k1=fhandler(x(i),y(i));;
    k2=fhandler(x(i)+h/2,y(i)+k1*h/2);
    k3=fhandler(x(i)+h/2,y(i)+k2*h/2);
    k4=fhandler(x(i)+h,y(i)+k3*h);
    y(i+1)=y(i)+h*((k1/6)+(2*k2/6)+(2*k3/6)+(k4/6));
  endfor
  plot(x,y)
endfunction
