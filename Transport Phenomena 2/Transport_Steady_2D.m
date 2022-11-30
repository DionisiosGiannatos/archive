#Transport Phenomena 2: Heat and Mass Transfer Christmas Project
#Heavily modified by me, original by the lab teachers.

% This code solves the heat equation in a 2D plate,
% using the Finite Difference Method.
% 1) Steady state heat transfer
% 2) No heat generation
% 3) Constant thermal conductivity
% 4) Fixed boundry conditions (Dirichlet Conditions)

clear;clc;hold off;close all; tic;

fprintf("Computing the Steady State of a 2D plate with two materials in contact \n\n")

fprintf("Starting the script: \n")
fprintf("Setting the initial parameters... \n")
% Material properties
k_a = 5; %Thermal conductivity 1
k_b = 1.5; %Thermal Conductivity 2

Lx = 4*10^(-2); % 1st plate dimension in x direction
Ly = 2*10^(-2); % 1st plate dimension in y direction
%-------------------------------------------------------
% Define the boundary conditions (Dirichlet Conditions)
T_top    = 25 ; % temperature on the upper side ( at y=Ly )
T_bottom = 25 ; % temperature on the lower side ( at y=0  )
T_left   = 100 ; % temperature on the left  side ( at x=0  )
           % temperature on the right side ( at x=Lx ) 
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
fprintf("\nCreating the linear system...\n")
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
for i=B
    k=k_a;
end

for i=middle
    k=(k_a+k_b)/2;
end

for i=C
    k=k_b;
end
% Plot the results
figure(2);plotContours(x,y,T);hold on;
plotGradient(x,y,T,dh,k)

hold off



