%{
This code works and it produces a graph but I'm not 100% sure 
what the graph includes. 
I think the x-axis is time and the y-axis is position. 
%}

%varying magnetic field
%no electric field
%ode23
clear;



x0 = 0; %starting position (m) analytic solution
y0 = 0; %starting position (m) analytic solution
z0 = 0;
vx = 1; %velocity in the x direction (m/s) analytic solution
vy = 0; %velocity in the y direction (m/s) ana solution
vz = 0;

KE_ana = sqrt(vx(1)^2 + vy(1)^2);

tspan = [0 50];
z0 = [vx,vy,vz,x0,y0,z0];
B0 = 1;
Bz = 1; %magnetic field in Tesla (T)
Bx = 0;
By = 0;
m = 1;%9.1038*10^(-31);% mass of electron (kg)
q = 1;%-1.602*10^(-19);% charge of electron; Coulombs

L = Bz*q/m; %omega for analytic solution ... not sure why I calculated this
T = (2*pi*m)/(abs(q)*Bz);% cyclotron period
vperp_o = sqrt(vx(1)^2 + vy(1)^2);%velocity perpendicular to magnetic field
Rc = m*vperp_o/(abs(q)*Bz);%radius of curvature

[t,y] = ode23( @(t,y) ode23magvar(t,y,B0,Bx,By,Bz,q,m) , tspan, z0);


plot(t,y(:,5))

function [ dzdt ] = ode23magvar( t,z,B0,Bx,By,Bz,q,m )
%UNTITLED5 for constant magnetic field 
%/Users/Angel/Documents/MATLAB/van allen radiation belt  project/ode23/ode23_mag_trial1.m
%/Users/Angel/Documents/MATLAB/van allen radiation belt  project/ode23/ode23_mag_trial1_figures.m
dzdt(1,1) = q/m*(z(2)*(1+z(4)/10)*B0-z(3)*By);%zzzz*1/(1+z(4)^2);
dzdt(2,1) = q/m*(z(3)*Bx-z(1)*(1+z(4)/10)*B0);
dzdt(3,1) = 0;
dzdt(4,1) = z(1);
dzdt(5,1) = z(2);
dzdt(6,1) = z(3);

% dzdt(1,1) = (z(2)*(1+z(4)/10)*B0);%zzzz*1/(1+z(4)^2);
% dzdt(2,1) = -z(1)*(1+z(4)/10)*B0;
% dzdt(3,1) = 0;
% dzdt(4,1) = z(1);
% dzdt(5,1) = z(2);
% dzdt(6,1) = z(3);

end
