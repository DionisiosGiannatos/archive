function output = least_square (x,y)
#An attempt at the creation of a general least-square script

#x: the x elements
#y: the y elements
#n: the amount of n elements

n_x=length(x);
n_y=length(y);

if n_x ~= n_y;
  error("x has a different amount of elements than y")
endif

n=n_x;

 SUM_x=0;
 SUM_y=0;
 xy=0;
 x_sqr=0;
 SSR=0;
 SSE=0;
 TSS=0;
 y_calc=0;
 
 for i=1:n;
   SUM_x=SUM_x + x(i);
   SUM_y=SUM_y + y(i);
   xy=xy+x(i)*y(i);
   x_sqr=x_sqr+x(i)^2;
 endfor
  x_ave=SUM_x./n;
  y_ave=SUM_y./n;
  
  xy_ave=x_ave*y_ave;
  x_sqr_ave=x_ave^2;
  
  a=(xy-n*xy_ave)/(x_sqr-n*x_sqr_ave);
  b=y_ave-a*x_ave;
  
  for i=1:n
    y_calc(i)=a*x(i)+b;
    SSR=SSR+(y_calc(i)-y_ave)^2;
    SSE=SSE+(y(i)-y_calc(i))^2;
    TSS=TSS+(y(i)-y_ave)^2;
  endfor
  
  R_sqrd=(1-(SSE/TSS));
  disp("Output: Alpha, Beta, R Squared");
  output=[a,b,R_sqrd]
endfunction
