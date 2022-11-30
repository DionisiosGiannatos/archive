function [g1, g2] = Van_Laar (x1_old, y1_old, T, P, x1, A1, B1, C1, A2, B2, C2)
  #A script for calculating the coefficients of the Van Laar equation
  # calculating g for known composition
  x2_old = 1-x1_old
  y2_old = 1-y1_old
  x2 = 1-x1
  P1s = Antoine (T, A1, B1, C1)
  P2s = Antoine (T, A2, B2, C2)
  g1_old = (y1_old*P)/(x1_old*P1s)
  g2_old = (y2_old*P)/(x2_old*P2s)

  #Van Laar Coeffs
  G = ((1+(x2_old*log(g2_old))/(x1_old*log(g1_old)))^2)*log(g1_old)
  D = ((1+(x1_old*log(g1_old))/(x2_old*log(g2_old)))^2)*log(g2_old)

  #G for unknown composition
  lng1 = G/(1+(G*x1)/(D*x2))^2
  lng2 = D/(1+(D*x2)/(G*x1))^2
  g1 = exp(lng1)
  g2 = exp(lng2)

endfunction