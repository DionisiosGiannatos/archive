#Transport Phenomena 2: Heat and Mass Transfer Christmas Project
#Heavily modified by me, original by the lab teachers.

% This code solves the heat equation in a 2D plate,
% using the Finite Difference Method.
% 1) Steady state heat transfer
% 2) No heat generation
% 3) Constant thermal conductivity
% 4) Fixed boundry conditions (Dirichlet Conditions)

clear;clc;hold off;close all; tic;
fprintf("Computing the Transient State of a 2D plate with two materials in contact \n\n")

fprintf("Starting the script: \n")
fprintf("Setting the initial parameters... \n")

% Material properties
k_a = 5; %Thermal conductivity 1
k_b = 1.5; %Thermal Conductivity 2
rho=7500; %Material density. Both materials have equal densities
cp=150; %Material Thermal Capacity, both materials have the same cp

transient_flag=false; %A flag used to determine whether we have reached the transient state or not. Will be used later.

a_a=k_a/(rho*cp); %Alpha factor of the A material
a_b=k_b/(rho*cp); %Alpha factor of the B material

Lx = 4*10^(-2); % 1st plate dimension in x direction
Ly = 2*10^(-2); % 1st plate dimension in y direction
%-------------------------------------------------------
% Define the boundary conditions (Dirichlet Conditions)
T_top    = 25 ; % temperature on the upper side ( at y=Ly )
T_bottom = 25; % temperature on the lower side ( at y=0  )
T_left   = 100 ; % temperature on the left  side ( at x=0  )
%-------------------------------------------------------


%--------------------------------------
%         Define the mesh            
%--------------------------------------
% Given dh (delta h)
dh=0.2*10^(-2); 
% Assumption: delta h = delta x = delta y 
% Define the (x,y) coordinates
x=0:dh:Lx; 
x(end)=Lx;
y=0:dh:Ly; 
y(end)=Ly;


nx=length(x);  % number of nodes in x direction
ny=length(y);  % number of nodes in y direction
n=nx*ny;       % total number of nodes 
% The dimension of the coefficient matrix A is n x n
% The dimension of the right hand side vector b is n x 1
%-------------------------------------
% Plot mesh
fprintf("\nPlotting the mesh... \n")
figure(1);plotMesh(x,y)
fprintf("Time elapsed in the plotting of the mesh:\n")
toc;
%--------------------------------------
% Mark nodes 
boundary=112:120;
right=(nx-1)*ny+1:(n-1);
left=1:ny;
top=ny:ny:n;
bottom=1:ny:nx*ny;

inner=1:n;
inner([left right top bottom boundary])=[];
%--------------------------------------


%--------------------------------------
%--------------------------------------
% Build the matrix A and the vector b %
%--------------------------------------
fprintf("\nCalculating the initial conditions... \n")
A=sparse(n,n); % Initialize the matrix A
b=zeros(n,1);  % Initialize the vector

for i=inner
    A(i,i-ny)=-1;
    A(i,i-1)=-1;
    A(i,i)=4;
    A(i,i+1)=-1;
    A(i,i+ny)=-1;
    
    b(i)=0;
end
% Impose Boundary conditions
for i=top
 %Dirichlet Condition: T_top=25 Celsius
    A(i,:)=0; A(i,i)=1;
    b(i)=T_top;
end

for i=bottom
 %Dirichlet Condition: T_top=25 Celsius
    A(i,:)=0; A(i,i)=1;
    b(i)=T_bottom;
end

for i=left
 %Dirichlet Condition: T_left=100 Celsius
    A(i,:)=0; A(i,i)=1;
    b(i)=T_left;
end
for i=right
 %Neumann Condition: Insulated Right Boundary
  A(i,i)=4;
  A(i,i-ny)=-2;
  A(i,i-1)=-1;
  A(i,i+1)=-1;
  
  b(i)=0;
end

  for i=boundary
    %Boundary between two materials with different conductivities 
    A(i,i)=2*(k_a+k_b);
    A(i,i-ny)=-k_a;
    A(i,i+ny)=-k_b;
    A(i,i-1)=-((k_a/2)+(k_b/2));
    A(i,i+1)=-((k_a/2)+(k_b/2));
    
    b(i)=0;
  end

%--------------------------------------


% Solve the linear system
T=A\b;
disp("Time for the solution of the  linear system:")
toc;

B=1:((nx-1)/2)*ny;
middle=((nx-1)/2)*ny+1:((nx+1)/2)*ny;
C=((nx+1)/2)*ny+1:n;
k=1:n;
for i=B
    k(i)=k_a;
endfor

for i=middle
    k(i)=(k_a+k_b)/2;
endfor

for i=C
    k(i)=k_b;
endfor
a=zeros(1,n);
for i=B
    a(i)=a_a;
endfor
for i=middle
  a(i)=k(i)/(rho*cp);
endfor
for i=C
    a(i)=a_b;
endfor
%--------- TRANSIENT STATE ---------
%Computation of the transient state using the Forward Euler method. Transient state occurs due to
%the increase of the left wall temperature at 130 degrees celsius (from 100 degrees).
%Initial condition is the solution of the previous problem.

%Step 1: Selection of timestep. The timestep was chosen arbitratily to be equal with 1 second, as to 
%not increase computational cost. In case of unstable behavior, the timestep has to be lowered at the
%cost of increased computational cost
fprintf("\nCalculation of the Transient State \n\n")
fprintf("Setting the initial parameters... \n")

dt=0.05;
t_final=50;
t=0; %initial time.
%timestep was chosen due to the accuracy it offered. any timestep bigger than approximately 0.45 produces
%wrong results.

%Step 2: Redefinition of the boundary conditions. The only condition that changes is the
%left wall temperature. Right wall remains insulated, top and bottom remain at 25 degrees celsius.
T_left_2   = 130 ; %Dirichlet condition at left wall.

%Step 3: Construction of the old temperature array. Mesh remains the same, only change is_absolute_filename
%Only change here is in the T_old vector, where the left-wall temperature is set to 130 degrees.
Told=T;
for i=left
 %Setting the left wall's temperature as 130 degrees celsius in order to obtain correct results.
 %Otherwise, the first loop will ignore the change in boundary conditions.
    Told(i)=T_left_2;
end

%Step 4: Actual calculation of the transient state.

%Initialization of the timeloop. This is where all calculations will be done.
fprintf("\nCalculating the transient state... \n")
for m=1:(t_final/dt)
  t=t+dt;
  for i=top
   %Dirichlet Condition: T_top=25 Celsius
    T(i)=T_top;
  end

  for i=bottom
   %Dirichlet Condition: T_top=25 Celsius
    T(i)=T_bottom;
  end

  for i=left
   %Dirichlet Condition: T_left=130 Celsius
    T(i)=T_left_2;
  end
  
  for i=boundary
    %Boundary between two materials with different conductivities 
   T(i)=Told(i)*(1-((2*(k_a+k_b)*dt)/(rho*cp*(dh^2))))+Told(i-ny)*(k_a*dt/(rho*cp*(dh^2)))+Told(i+ny)*(k_b*dt/(rho*cp*(dh^2)))+Told(i-1)*((k_a+k_b)/2)*(dt/(rho*cp*(dh^2)))+Told(i+1)*((k_a+k_b)/2)*(dt/(rho*cp*(dh^2)));
  end
  
  for i=inner
    %Calculation of temperature at an internal point. Equation originated from energy-balance.
    %In order to accurately represent the transfer on the actual material, the a constant on the west and east 
    %has been set to be equal with the a of the point. This produces more accurate results 
    %near the boundary rather than the hypothesis that the a constant of the boundary is equal to (ka+kb)/2
    T(i)= Told(i)*(1-((4*a(i)*dt)/(dh^2))) + Told(i-ny)*((a(i)*dt)/(dh^2)) + Told(i+ny)*((a(i)*dt)/(dh^2)) + Told(i-1)*((a(i-1)*dt)/(dh^2)) + Told(i+1)*((a(i+1)*dt)/(dh^2)) ;
  end
  
  for i=right
   %Neumann Condition: Insulated Right Boundary
    T(i)= Told(i)*(1-(4*a(i)*dt)/(dh^2)) + Told(i-ny)*((2*a(i-ny)*dt)/(dh^2)) + Told(i-1)*((a(i-1)*dt)/(dh^2)) + Told(i+1)*((a(i+1)*dt)/(dh^2)) ;
  end

  
  c1=norm(T); %calculation of the norms of the new and the old Temperature vectors. This will be used
  c2=norm(Told); %to determine when we have reached the transient state.
  criterion=abs((c1-c2)/c1); %The criterion of whether we have reached transient state or not.
  
  %A note on the steady-state calculation. The script works by detecting whether the changes to the 
  %vector's norm are significant or not. This means that a few elements might change by 1-2 degrees
  %Celsius, but because the impact to the norm is so low, the program believes that we have reached
  %a steady state.
  
  if (criterion<2*10^(-5) && transient_flag==false)
    t_steady=t; %captures the moment where the criterion became true.
    transient_flag=true; %a marker that we have reached the steady state.
    fprintf("Steady state reached approximately at t= %f seconds\n", t_steady)
  endif
  
  Told=T;
  
  if (m==200)
    disp("Temperature at boundary for 10 seconds. Temperature starts from node 112 to node 120");
    T(112:120)
    fprintf("Script paused. Press any key to continue...\n\n")
    pause();
   endif
 
  if (m==400)
    disp("Temperature at boundary for 20 seconds. Temperature starts from node 112 to node 120");
    T(112:120) 
    fprintf("Script paused. Press any key to continue...\n\n")
    pause();
  endif
  
  if (m==1000)
    disp("Temperature at boundary for 50 seconds. Temperature starts from node 112 to node 120");
    T(112:120) 
  endif

endfor 
fprintf("\nThe script has finished executing \n")
disp("Total Runtime:")
toc;
hold off



