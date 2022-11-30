#A simple Newton-Raphson Calculator. f must be an anonymous function.
function output = NR (f,xinit)
h=sqrt(eps)*abs(xinit);
fder=@(x,h) (f(x+h)-f(x-h))/(2*h); #Central Difference to calculate derivative.
x=xinit-(f(xinit)/fder(xinit,h));
xnew=x;
xold=xinit;
itr=0;
chk=abs(xnew/(xold)-1);
while (chk>10^(-8) && itr<200)
  itr=itr+1; #iteration counter.
  xold=x;
  h=sqrt(eps)*abs(x);
  x=x-(f(x)/fder(x,h));
  xnew=x;
  chk=abs(xnew/(xold)-1);
endwhile
if (itr==200)
  disp("Last output:")
  disp(xnew);
  error("No solution could be found. Try different initial guess, or check function.")
endif
disp("Output type: Solution, Iterations")
output=[xnew,itr]
endfunction