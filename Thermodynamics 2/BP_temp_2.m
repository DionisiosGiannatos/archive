function [Tf, y1, y2] = BP_temp_2(P, x1, A1, B1, C1, A2, B2, C2, g1, g2)
  x2=1-x1;
  # Boiling Point Calculation Tb
  disp("Boiling Points ")
  Tb1 = (B1/(A1-log(P))) -C1
  Tb2 = (B2/(A2-log(P))) -C2
  disp("Initial Temperature")
  T0 = x1*Tb1+x2*Tb2

  # Pressure Calculation for said Temperature
  P1s = Antoine (T0, A1, B1, C1)
  P2s = Antoine (T0, A2, B2, C2)
  P_temp = x1*g1*P1s + x2*g2*P2s
  if P_temp > P
    disp("Replace with the Maximum Temperature")
  else
    disp("Replace with the Minimum Temperature")
  endif
  Tbnew = input("Give the 2 Required Temperatures (Vector Form)");

  for i=1:20
    T = x1*Tbnew(1)+x2*Tbnew(2)
    P1s = Antoine (T, A1, B1, C1)
    P2s = Antoine (T, A2, B2, C2)
    P_temp = x1*g1*P1s + x2*g2*P2s
    if P_temp > P
      disp("Pressure calculated higher: Give pair of lowest temps")
    else
      disp("Pressure Calculated lower: give pair of highest temps")
    endif
    cont = input("Continue Loop? (0=no): ");
    if cont == 0
      printf("Loop stopped after %d iterations \n", i)
      disp("Pressure Error:")
      P_err = (abs(P_temp-P)/P)*100
      break
    endif
    Tbnew = input("Give me the 2 required Temperatures ");
  endfor

  Tf = T;
  y1 = (x1*g1*P1s)/P;
  y2 = (x2*g2*P2s)/P;
endfunction