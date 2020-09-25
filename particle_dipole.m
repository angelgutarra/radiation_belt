%A postively charged particles earth model as a magnetic dipole
%no electric field
%ode23
%The maximum number of times the process can go is 460 cyclotron periods. 
%Continuing the process beyond 460 periods causes the particle to be
%ejected from it's orbit. 
%At first I thought it was because the particle kept building speed until
%it was launched out orbit but after graphing the kinetic energy vs time it
%was clear that the particle was losing energy. 
%The pitch angle from initial velocities is 63.4 degrees. 
%The pitch angle calculated using the magntic field and the strongest point
%of the magnetic field yeilds 28 degrees. 
%The pitch angle has to be greater than 28 degrees for the particle to not
%be lost to the atmosphere of the Earth. 
clear;

B0 = 6.5*10^(-5); %earths' magnetic field at the surface of the Earth (T)
Re = 6.371*10^6; % radius of the Earth (m)

m = 1;%1.67*10^(-10);%9.1038*10^(-31);% mass of electron (kg)
q = 1;%1.602*10^(-19);% charge of electron; Coulombs

x0 = 2*Re; %starting position (m)
y0 = 0; %starting position (m) 
z0 = 0;
vx = 10; %velocity in the x direction (m/s) 
vy = 0; %velocity in the y direction (m/s) 
vz = 5;%velocity in the z direction (m/s) 

pitch_angle = atan(sqrt(vx(1)^2+vy(1)^2)/vz)*180/pi;%pitch angle 

% location of particle in spherical coordinates
r = sqrt(x0^2+y0^2+z0^2);
TEA = acos(z0/r);
PHI = atan(y0/(x0));

% magnetic field in spherical coordinates at the location of the particle   
Br = 2*B0*(Re/r)^3*cos(TEA);
BTEA = B0*(Re/r)^3*sin(TEA);
    
% magnteic field in rectangular coordinates
B = sqrt(Br^2 + BTEA^2);
Bz = Br*cos(TEA)-BTEA*sin(TEA);
Bxt = Br*sin(TEA)+BTEA*cos(TEA);
Bx = Bxt*cos(PHI);
By = Bxt*sin(PHI);

T = (2*pi*m)/(abs(q)*B);% cyclotron period
vperp_o = sqrt(vx(1)^2 + vy(1)^2 );%velocity perpendicular to magnetic field
Rc = m*vperp_o/(abs(q)*B);%radius of curvature

KE = sqrt(vx(1)^2 + vy(1)^2 + vz(1)^2);

tspan = [0 450*T];
z0 = [vx,vy,vz,x0,y0,z0];
[t,y] = ode23( @(t,y) ode23dipole(t,y,Bx,By,Bz,q,m,B0,Re) , tspan, z0);
plot3(y(:,4)/Re,y(:,5)/Re,y(:,6)/Re)
grid on
title('A Positvely Charged Particle''s Trajectory in Earths Magntictic Field Modeled as a Dipole')
xlabel('x postion (x/Re)')
ylabel('y postion (y/Re)')
zlabel('z position (z/Re)')


for i = 1:length(t)
    %Kinetic energy
    KE(i) = 1/2*m*(y(i,1)^2+y(i,2)^2+y(i,3)^2);
    
    
    %%strength of magnetic field at the location of the particle

    % location of particle in spherical coordinates
    r(i) = sqrt(y(i,4)^2+y(i,5)^2+y(i,6)^2);
    TEA(i) = acos(y(i,6)/r(i));
    PHI(i) = atan(y(i,5)/(y(i,4)));

    % magnetic field in spherical coordinates at the location of the particle   
    Br(i) = 2*B0*(Re/r(i))^3*cos(TEA(i));
    BTEA(i) = B0*(Re/r(i))^3*sin(TEA(i));
    
    % magnteic field in rectangular coordinates
    B(i) = sqrt(Br(i)^2 + BTEA(i)^2);
end

%pitch_angle according to computational method
pitch_a_comp = asin(sqrt(B(1)/max(B)))*180/pi;
