function [V, L, x, y] = flash (K, z, F)
  #One of the first applications of the Newton-Raphson method, used for the calculation of a Flash Evaporation

  #K,z,F are vectors
  # Variable Initialization
  c = length(z);
  a=0.5;
  for i=1:c
    Q(i) = (z(i)*(K(i)-1))/(1+a*(K(i)-1));
    Qder(i) = -(Q(i)^2)/z(i);
  endfor
  disp("Q: ")
  disp(Q)
  disp("Derivative of Q: ")
  disp(Qder)
  sumQ = sum(Q)
  sumQder = sum(Qder)

  # Newton-Raphson
  for n=1:3
    disp("fraction of Q and Derivative of Q")
    disp(sumQ/sumQder)
    a = a - (sumQ/sumQder)
    for i=1:c
      Q(i) = (z(i)*(K(i)-1))/(1+a*(K(i)-1));
      Qder(i) = -(Q(i)^2)/z(i);
    endfor
    disp("Q: ")
    disp(Q)
    disp("Derivative of Q: ")
    disp(Qder)
    sumQ = sum(Q)
    sumQder = sum(Qder)
  endfor 

  # Final Result:
  disp("The final value of a is:")
  disp(a)
  for i=1:c
    printf("Component: %d \n", i)
    x(i)=z(i)/(1+a*(K(i)-1))
    y(i)=(K(i)*z(i))/(1+a*(K(i)-1))
  endfor
  disp("Amount of Vapor:")
  V=a*F
  disp("Amount of Liquid:")
  L=F-V

endfunction 