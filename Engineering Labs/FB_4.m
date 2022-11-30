function y = FB_4(x) 
% Constants: These remain constant throughout the experiment
  g=9.81; %gravity constant m/s^2
  rho_s= 2900; %density of fluidized material, kg/m^3
  rho=1.213; %density of fluid, kg/m^3
  d_eq=0.000789157; %equivalent diameter of material, m
  dyn_vis=0.0000181; %dynamic viscosity, Pa*s
  l_star=0.039085; %initial bed length, m
  e_star=0.640118; %initial porosity

% Variables: These change depending on the air feed
  Dp=540; %pressure drop of column, Pascal
  u_phi= 2.72859118; %apparent velocity of fluid, m/s
  
% Actual Solution by octave's fzero. x(1) is the equivalent porosity, x(2) is the equivalent length
  y=zeros(2,1);
  y(1)=l_star*((1-e_star)/(1-x(1)))-x(2);
  y(2)=(x(2)/d_eq)*((1-x(1))/(x(1)^(3)))*(1.75*(u_phi^2)+((150*(1-x(1))*dyn_vis*u_phi)/(rho*d_eq)))-(Dp/rho);
endfunction