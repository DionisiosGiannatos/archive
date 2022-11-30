#A forward euler script, created for practice
#requires a function dy/dt=f(y,t), and an initial value.
#vinit and vfinal represent the initial and final values.
#npoints refers to the number of points
function retval = f_euler (f,yinit,vinit=0,vfinal=10,npoints=1000)
  a=linspace(vinit,vfinal,npoints);
  h=a(2)-a(1);
  y(1)=yinit;
  plot(a(1),y(1))
  for i=1:(npoints-1);
    y(i+1)=yinit+f(y,a(i))*h
    plot(a(i+1),y(i+1))
  endfor
endfunction
