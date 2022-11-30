#Made for practice for the Numerical Methods class
#ex_finder stands for extrema finder.
#Extrema is the terminology utilized to describe Maximum or Minimum of a function.
#Inputs is an anonymous function f and an initial guess
#Central Differences used to find 1st and 2nd order derivatives
#Funny fact about polyonomials: trying x^2 doesn't give a root cause it actually finds it 1st try.
#Check however tells program to keep going cause it results in (0/1)-1=1
function output = ex_finder (f, xinit)
h=sqrt(eps)*abs(xinit);
f_1der=@(x,h) (f(x+h)-f(x-h))/(2*h);
f_2der=@(x,h) (f(x+h)-2*f(x)+f(x-h))/(h^2);
x=xinit-(f_1der(xinit,h)/f_2der(xinit,h));
xnew=x;
xold=xinit;
itr=0;
chk=abs(xnew/(xold)-1);
while (chk>10^(-8) && itr<200)
  itr=itr+1; 
  xold=x;
  h=sqrt(eps)*abs(x);
  x=xold-(f_1der(xold,h)/f_2der(xold,h));
  xnew=x;
  chk=abs(xnew/(xold)-1);
endwhile
if (itr==200)
  disp("Last output:")
  disp(xnew);
  error("No solution could be found. Try different initial guess, or check function.")
endif
if (f_2der(xnew,h)<0)
  disp("Maximum Located")
elseif (f_2der(xnew,h)>0)
  disp("Minimum Located")
endif
disp("Output type: Solution, Iterations")
output=[xnew,itr]
endfunction
