function output = e_calc (k1,k2)
#I am not sure what this script does, or even which class it was for. I am not even sure if its functional
e1=0;
e2=0;
iter=0;
diff=10;
e1_correct=0;
e2_correct=0;
root_found=logical(false);

#equation solving k2 as e1, calculating e1 through a preexisting e2: (-3*k2)e1^(2)+e1(3*k2+(-2*k2-1)*e2) +((k2-1)*e2^(2)-k2*e2)
#equation solving k1 as e2, calculating e2 through a preexisting e1: e2^(4)-8*e1*e2^(3)+e2^(2)*18*e1^(2)+e2(-20.4*e1^(3)-20.4*e1^(2)+35.7*e1-10.2)+(10.2-45.9*e1-47.4*e1^(4)-56.1*e1^(2))
while root_found == false && e1<1;
 e1=e1+0.0001;
 iter=iter+1;
 e2_pre=roots([1,-8*e1,18*e1^(2),(-20.4*e1^(3)-20.4*e1^(2)+35.7*e1-10.2),(10.2-45.9*e1-47.4*e1^(4)-56.1*e1^(2))]);
 e2_pre=e2_pre(e2_pre==real(e2_pre));
 n=rows(e2_pre);
 for i=1:n
   if (e2_pre(i,1) < 1) && (e2_pre(i,1)>0)
     e2=e2_pre(i,1);
   else
     continue 
   endif
   e1_pre=roots([(-3*k2),(3*k2+(-2*k2-1)*e2),((k2-1)*e2^(2)-k2*e2)]);
   e1_pre=e1_pre(e1_pre==real(e1_pre));
   k=rows(e1_pre);
   if k=0;
     continue
   endif
   for j=1:k
     if (e1_pre(j,1) < 1) && (e1_pre(j,1) > 0);
       e1_calc=e1_pre(j,1);
       diff=abs((e1/e1_calc)-1);
     else
       continue
     endif
     if diff<10^(-4);
       root_found=logical(true);
       e2_correct=e2;
       e1_correct=e1_calc;
     endif
    endfor
   endfor
endwhile 
if e1=1;
  error("Root not calculated")
endif
output=[e1_correct,e2_correct,iter];
endfunction
