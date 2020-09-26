%constant mag field
%figures for ode23_mag_trial1.m 
%/Users/Angel/Documents/MATLAB/van allen radiation belt  project/ode23/ode23_mag_trial1.m
%/Users/Angel/Documents/MATLAB/van allen radiation belt  project/ode23/ode23fn.m
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
B0 = 5;
Bz = 1; %magnetic field in Tesla (T)
Bx = 1;
By = 1;
m = 1;%9.1038*10^(-31);% mass of electron (kg)
q = 1;%-1.602*10^(-19);% charge of electron; Coulombs

L = Bz*q/m; %omega for analytic solution ... not sure why I calculated this
T = (2*pi*m)/(abs(q)*Bz);% cyclotron period
vperp_o = sqrt(vx(1)^2 + vy(1)^2);%velocity perpendicular to magnetic field
Rc = m*vperp_o/(abs(q)*Bz);%radius of curvature

[t,y] = ode23( @(t,y) ode23fn(t,y,Bx,By,Bz,q,m) , tspan, z0);


i = 1;

%postions versus time x,y,z
figure (i)
plot(t,y(:,4))
grid on
title ('constant magnetic field in the z+ direction')
xlabel('time (s)')
ylabel('x position (m)')
i = i+1;

figure (i); clf
plot(t,y(:,5))
grid on
title ('constant magnetic field in the z+ direction')
xlabel('time (s)')
ylabel('y position (m)')
i = i+1;

figure (i); clf
plot(t,y(:,6))
grid on
title ('constant magnetic field in the z+ direction')
xlabel('time (s)')
ylabel('z position (m)')
i = i+1;

%velocity vs time
figure (i); clf
plot(t,y(:,1))
grid on
title ('constant magnetic field in the z+ direction')
xlabel('time (s)')
ylabel('x velocity (m/s)')
i = i+1;

figure (i); clf
plot(t,y(:,2))
grid on
title ('constant magnetic field in the z+ direction')
xlabel('time (s)')
ylabel('y velocity (m/s)')
i = i+1;

figure (i); clf
plot(t,y(:,3))
grid on
title ('constant magnetic field in the z+ direction')
xlabel('time (s)')
ylabel('z velocity (m/s)')
i = i+1;


%postion
figure(i);clf
plot(y(:,4),y(:,5))
grid on
title('constant magnetic field in the z+ direction')
xlabel('x postion (m)')
ylabel('y postion (m)')
i = i+1;

figure(i);clf
plot(y(:,5),y(:,6))
grid on
title('constant magnetic field in the z+ direction')
xlabel('y postion (m)')
ylabel('z postion (m)')
i = i+1;

figure(i);clf
plot(y(:,4),y(:,6))
grid on
title('constant magnetic field in the z+ direction')
xlabel('x postion (m)')
ylabel('z postion (m)')
i = i+1;


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
