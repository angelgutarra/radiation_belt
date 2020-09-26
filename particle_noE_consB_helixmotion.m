%ode 23
%constant mag field
%no elec field
%/Users/Angel/Documents/MATLAB/van allen radiation belt  project/ode23/ode23fn.m
%/Users/Angel/Documents/MATLAB/van allen radiation belt  project/ode23/ode23_mag_trial1_figures.m

Bz = 1; %magnetic field in Tesla (T)
Bx = 0;
By = 0;
m = 1;%9.1038*10^(-31);% mass of electron (kg)
q = 1;%-1.602*10^(-19);% charge of electron; Coulombs
L = Bz*q/m; %omega for analytic solution

x0 = 0; %starting position (m) analytic solution
y0 = 0; %starting position (m) analytic solution
z0 = 0;
vx = 1; %velocity in the x direction (m/s) analytic solution
vy = 0; %velocity in the y direction (m/s) ana solution
vz = 1;

T = (2*pi*m)/(abs(q)*Bz);% cyclotron period
vperp_o = sqrt(vx(1)^2 + vy(1)^2);%velocity perpendicular to magnetic field
Rc = m*vperp_o/(abs(q)*Bz);%radius of curvature

KEana = sqrt(vx(1)^2 + vy(1)^2);

tspan = (0:0.1:10);
z0 = [vx,vy,vz,x0,y0,z0];
Bz=1;
[t,y] = ode23( @(t,y) ode23fn(t,y,Bx,By,Bz,q,m) , tspan, z0);


i = 1;
Re = 1;
%postions versus time x,y,z
figure (i)
plot(t/T,y(:,4)/Re)
grid on
title ('A Postively Charged Particles Postition in x-direction vs. time')
xlabel('time (t/T)')
ylabel('x position (x/Re)')
i = i+1;

figure (i); clf
plot(t/T,y(:,5)/Re)
grid on
title ('A Postively Charged Particles Position in the y-direction vs. Time')
xlabel('time (t/T)')
ylabel('y position (y/Re)')
i = i+1;

figure (i); clf
plot(t/T,y(:,6)/Re)
grid on
title ('A Postively Charged Particles Position in the z-direction vs. Time')
xlabel('time (t/T)')
ylabel('z position (z/Re)')
i = i+1;

%velocity vs time
figure (i); clf
plot(t/T,y(:,1))
grid on
title ('A Postively Charged Particles Velocity in the x-direction vs. Time')
xlabel('time (t/T)')
ylabel('x velocity (m/s)')
i = i+1;

figure (i); clf
plot(t/T,y(:,2))
grid on
title ('A Postively Charged Particles Velocity in the y-direction vs. Time')
xlabel('time (t/T)')
ylabel('y velocity (m/s)')
i = i+1;

figure (i); clf
plot(t/T,y(:,3))
grid on
title ('A Postively Charged Particles Velocity in the z-direction vs. Time')
xlabel('time (t/T)')
ylabel('z velocity (m/s)')
i = i+1;


%postion relationship
figure(i);clf
plot(y(:,4)/Re,y(:,5)/Re)
grid on
title('A Postively Charged Particles Position in the x-direction vs y-direction')
xlabel('x postion (x/Re)')
ylabel('y postion (y/Re)')
i = i+1;

figure(i);clf
plot(y(:,5)/Re,y(:,6)/Re)
grid on
title('A Postively Charged Particles Position in the y-direction vs z-direction')
xlabel('y postion (y/Re)')
ylabel('z postion (z/Re)')
i = i+1;

figure(i);clf
plot(y(:,4)/Re,y(:,6)/Re)
grid on
title('A Postively Charged Particles Position in the x-direction vs z-direction')
xlabel('x postion (x/Re)')
ylabel('z postion (z/Re)')
i = i+1;

%3D plot
figure(i);clf
plot3(y(:,4)/Re,y(:,5)/Re,y(:,6)/Re)
grid on
title('Earth as a dipole')
xlabel('x postion (x/Re)')
ylabel('y postion (y/Re)')
zlabel('z position (z/Re)')
i = i+1;

%Kinetic energy 
% this graph is not producing the right results
% What I was probably trying to do was graph 
% the KE from interpolated velocities in one plot and 
% then compare it to the KE I obtained from an analytic 
% equation. 
figure(i);clf
plot(t/T,KEana/KEana(1))
grid on
title('Kinetic Energy of the Particle vs. Time')
xlabel('time (t/T)')
ylabel('KE/KE(1)')

function [ dzdt ] = ode23fn( t,z,Bx,By,Bz,q,m )
%UNTITLED5 for constant magnetic field 
%/Users/Angel/Documents/MATLAB/van allen radiation belt  project/ode23/ode23_mag_trial1.m
%/Users/Angel/Documents/MATLAB/van allen radiation belt  project/ode23/ode23_mag_trial1_figures.m
dzdt(1,1) = q/m*(z(2)*Bz-z(3)*By);%zzzz*1/(1+z(4)^2);
dzdt(2,1) = q/m*(z(3)*Bx-z(1)*Bz);
dzdt(3,1) = 0;
dzdt(4,1) = z(1);
dzdt(5,1) = z(2);
dzdt(6,1) = z(3);

end

