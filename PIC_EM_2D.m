%%% This code is used for 3-D particle in cell plasma simulation
clear;
close all;
%%% Constants will be defined as below
global miu0 ep0
miu0 = 1.25663706127e-6 % N⋅A−2
ep0 = 8.8541878188e-12 % F⋅m−1
me = 9.1093837e-31 % kg
mi = 1.67e-27 % kg
e0 = 1.602176634e-19 % C
kb = 1.380649e-23 %J⋅K−1

%User set constant
T0_e = 1e1 %K
T0_i = 1e1 %K
ne = 1e14 %m-3
global grid_size
grid_size = [10,20] %number of grid in each direction
total_step = 100000 %total step of the simulation
global particle_number_i particle_number_e
particle_number_setted = 1e3 % Setted number of particle used in this simulation
particle_number_i = particle_number_setted %Number of ion used in this simulation
particle_number_e = particle_number_setted %Number of electron used in this simulation
v_drift = 0
max_injection = 10

%Calculated constant
global debye_length
debye_length = sqrt(ep0*kb*T0_e/(ne*e0^2)) %m
plasma_frequency = sqrt(ne*e0^2/(me*ep0)) %rad/s
global t_unit x_total y_total
t_unit = 0.1 / plasma_frequency %s time unit of the simutaion
x_total = grid_size(1) * debye_length %Physical x length simulated
y_total = grid_size(2) * debye_length %Physical x length simulated
total_number = x_total*y_total*ne %Total electron number
global number_represent
number_represent = total_number / particle_number_setted %Number of electrons a particle presents

e_p = -e0 * number_represent %Charge of a electron
m_p_e = me * number_represent %Mass of a electron
e_i = e0 * number_represent %Charge of a ion
m_p_i = mi * number_represent %Mass of a ion

v0_e = sqrt(kb*T0_e*2/(me*number_represent))
v0_i = sqrt(kb*T0_i*2/(mi*number_represent))

global d
d = debye_length^2/(2*t_unit)%E correction fector, based on
%Barry Marder,
%A method for incorporating Gauss' law into electromagnetic PIC codes,
%Journal of Computational Physics,
%Volume 68, Issue 1,
%1987,
%Pages 48-55,
%ISSN 0021-9991,
%https://doi.org/10.1016/0021-9991(87)90043-X.
%(https://www.sciencedirect.com/science/article/pii/002199918790043X)


% E solver
function E = E_solver(E0,B,J)
%%%
% This is the solver that solves E in next time step using maxwll equation
%
% Input
% E0 is a structure contains E0.x and E0.y which are gird_size(1) *
% gird_size(2) 2D arrays contains x components and y components of the E
% field
% B is a grid_size(1) * gird_size(2) 2D array contains z component of
% magnetic field
% J is a structure contains J.x and J.y which are gird_size(1) *
% grid_size(2) 2D arrays contains x and y components of current density J
%
% Output
% E is a structure contains E.x and E.y which are gird_size(1) *
% gird_size(2) 2D arrays contains x components and y components of the E
% field in next time step
%%%

global miu0 ep0 debye_length t_unit

[dx,dy] = gradient(B,debye_length);
E.x = E0.x + t_unit * (- J.x/(ep0));
E.y = E0.y + t_unit * (- J.y/(ep0));

end

% B solver
function B = B_solver(B0,E)
%%%
% This is the solver that solves B in next time step using maxwll equation
%
% Input
% B0 is a grid_size(1) * gird_size(2) 2D array contains z component of
% magnetic field
% E is a structure contains E0.x and E0.y which are gird_size(1) *
% gird_size(2) 2D arrays contains x components and y components of the E
% field
%
% Output
% B0 is a grid_size(1) * gird_size(2) 2D array contains z component of
% magnetic field in next time step
%%%

global t_unit debye_length
[xdx,xdy] = gradient(E.x,debye_length);
[ydx,ydy] = gradient(E.y,debye_length);

B = B0 - t_unit * (ydx-xdy);
end

% Weighting
function [x,y,w] = weighting(x0,y0)
%%%
% This function is used to calculate weighting for a particle located at
% (x0,y0)
% x0,y0 are x and y coordinates of the particle
% [x,y] are a number / a array of coordinates that particle contributes
% w is a number / array of weighing each cooridnate wights
%%%

global debye_length grid_size

xi = floor(x0/debye_length)+1;
dx = (x0/debye_length)+1 - xi;
yi = floor(y0/debye_length)+1;
dy = (y0/debye_length)+1 - yi;

x = [xi+1,xi,xi,xi+1];
y = [yi+1,yi+1,yi,yi];

w = [dx*dy,(1-dx)*dy,(1-dx)*(1-dy),dx*(1-dy)];

 %x = floor(x0/debye_length)+1;
 %y = floor(y0/debye_length)+1;
 %w = 1;

for i = 1:length(x)
    if x(i) > grid_size(1)
        x(i) = grid_size(1);
    end
    if y(i) > grid_size(2)
        y(i) = grid_size(2);
    end
    if x(i) < 1
        x(i) = 1;
    end
    if y(i) < 1
        y(i) = 1;
    end
end
w = w/sum(w);
end

% Particle motion
function [particle,velocity] = particle_motion(n,p0,v0,e,m,E,B)
%%%
% This function is used to calculate particle motion to next step
%
% Input
% n is the particle number
% p0 is a particle_number# * 2 array contains x and y coordinates of the
% particles
% v0 is a particle_number# * 2 array contains x and y components of the
% particles
% e is the charge of the gianent particle
% m is the mass of the gianent patrticle
% E is the a structure contains E.x and E.y that are gird_size(1)
% * gird_size(2) contains x and y components of the electronic field
% B is a gird_size(1) * gird_size(2) array contains z componenets of the
% magnetic field
%
% Output
% particle is a particle_number# * 2 array contains x and y coordinates
% of the particles in next time step
% velocity is a particle_number# * 2 array contains x and y components of
% the particles in next time step
%%%

global t_unit E_ext B_ext
particle = zeros(n,2);
velocity = zeros(n,2);
for i = 1:n
    B_local = 0;
    E_local.x = 0;
    E_local.y = 0;
    [x,y,w] = weighting(p0(i,1),p0(i,2));
    for j = 1:length(x)
        B_local = B_local + B(x(j),y(j)) * w(j) + B_ext(x(j),y(j)) * w(j);
        E_local.x = E_local.x + E.x(x(j),y(j)) * w(j) + E_ext.x(x(j),y(j)) * w(j);
        E_local.y = E_local.y + E.y(x(j),y(j)) * w(j) + E_ext.y(x(j),y(j)) * w(j);
    end
    % a.x = (e/m) * E_local.x + v0(i,2)*B_local;
    % a.y = (e/m) * E_local.y - v0(i,1)*B_local;
    %
    % velocity(i,1) = v0(i,1) + a.x * t_unit;
    % velocity(i,2) = v0(i,2) + a.y * t_unit;
    %
    % particle(i,1) = p0(i,1) + 0.5*(velocity(i,1)+v0(i,1)) * t_unit;
    % particle(i,2) = p0(i,2) + 0.5*(velocity(i,2)+v0(i,2)) * t_unit;

    %28 Jan 2025 : Here we try Boris push method
    vx_n = v0(i,1) + 0.5 * (e/m) * E_local.x * t_unit;
    vy_n = v0(i,2) + 0.5 * (e/m) * E_local.y * t_unit;

    vx_p = ((1-(t_unit*e*B_local/m)^2)*vx_n + vy_n*2*t_unit*e*B_local/m)/(1+(t_unit*e*B_local/m)^2);
    vy_p = ((1-(t_unit*e*B_local/m)^2)*vy_n - vx_n*2*t_unit*e*B_local/m)/(1+(t_unit*e*B_local/m)^2);

    velocity(i,1) = vx_p + 0.5 * (e/m) * E_local.x * t_unit;
    velocity(i,2) = vy_p + 0.5 * (e/m) * E_local.y * t_unit;

    particle(i,1) = p0(i,1) + 0.5*(velocity(i,1)+v0(i,1)) * t_unit;
    particle(i,2) = p0(i,2) + 0.5*(velocity(i,2)+v0(i,2)) * t_unit;
end
end

% Boundary conditions
function [n,particle,velocity] = boundary(n0,p0,v0)
%%%
% This function is used to determin boundary condition of the particles,
% i.e. either reflect, move it somewhere or kill it
%
% Input
% n0 is the number of particle
% p0 is a particle_number# * 2 array contains x and y coordinates of the
% particles
% v0 is a particle_number# * 2 array contains x and y components of the
% particles
%
% Output
% n is the new particle number
% particle is a n# * 2 array contains x and y coordinates
% of the particles in next time step
% velocity is a n# * 2 array contains x and y components of
% the particles in next time step
%%%

global x_total y_total
n = 0;
particle = zeros(n,2);
velocity = zeros(n,2);
for i = 1 : n0
    %reflect at x and y boundary
    if p0(i,1) > x_total
        p0(i,1) = x_total - (p0(i,1) - x_total);
        v0(i,1) = - v0(i,1);
    end
    if p0(i,1) < 0
        p0(i,1) = -p0(i,1);
        v0(i,1) = -v0(i,1);
    end
    if p0(i,2) > y_total
        p0(i,2) = y_total - (p0(i,2) - y_total);
        v0(i,2) = - v0(i,2);
    end
    if p0(i,2) < 0
        p0(i,2) = -p0(i,2);
        v0(i,2) = -v0(i,2);
    end

    if p0(i,1) >= 0 && p0(i,1) <= x_total &&  p0(i,2) >= 0 && p0(i,2) <= y_total
        particle(n+1,:) = p0(i,:);
        velocity(n+1,:) = v0(i,:);
        n = n+1;
    end
end
end

% Calculate J
function J = J_solver(n0,p0,v0,e)
%%%
% This function is used to solve current density J
%
% Input
% n0 is the number of particle
% p0 is a particle_number# * 2 array contains x and y coordinates of the
% particles
% v0 is a particle_number# * 2 array contains x and y components of the
% particles
% e is the charge of the giaent particle
%
% Output
% J is a structure contains J.x and J.y which are grid_size(1) *
% gird_size(2) array contains x and y component of J
%%%

global grid_size debye_length
J.x = zeros(grid_size(1),grid_size(2));
J.y = zeros(grid_size(1),grid_size(2));

for i = 1 : n0
    [x,y,w] = weighting(p0(i,1),p0(i,2));
    for j = 1:length(x)
        J.x(x(j),y(j)) = J.x(x(j),y(j)) + w(j) * e * v0(i,1) /debye_length^2;
        J.y(x(j),y(j)) = J.y(x(j),y(j)) + w(j) * e * v0(i,2) /debye_length^2;
    end
end
end

% Calculate total potential energy
function total_U = potential(n0_1,n0_2,p0_1,p0_2,q1,q2)
%%%
% This functiuon is used to calculate total potential enegry of the
% particles
%
% Input
% n0_1 is the number of particle
% p0_1 is a particle_number# * 2 array contains x and y coordinates of the
% particles
% n0_2 is the number of particle
% p0_2 is a particle_number# * 2 array contains x and y coordinates of the
% particles
% q1 is the charge in n0_1
% q2 is the chagre in n0_2
% Output
% total_U is the total potential energy 
%%%

global ep0 debye_length
total_U = 0;

for i = 1:n0_1
    for j = i+1:n0_2
        r = sqrt((p0_1(i,1)-p0_2(j,1))^2+(p0_1(i,2)-p0_2(j,2))^2);
        if r == 0 
            continue
        end
        total_U = total_U + (q1*q2)/(4*pi*ep0*r);
    end
end
total_U
end

% Density ref
function den = den_ref (E)
%%%
% This function is used to double check E value by calculate charge density
% with Gauss`s law
%
% Input
% E is a structure contains E0.x and E0.y which are gird_size(1) *
% gird_size(2) 2D arrays contains x components and y components of the E
% field
%
% Output
% den is a gird_size(1) * gird_size(2) 2D array contain charge density, SI
% unit
%%%

global ep0 debye_length

[xdx,~] = gradient(E.x,debye_length);
[~,ydy] = gradient(E.y,debye_length);
den =  ep0 * (xdx+ydy);
end

% Plot a figure
function den = density(n0,p0)
%%%
% This function is used to calculate number density
%
% Input
% n0 is the number of particles
% p0 is a n0*2 array contains position imformation of the particles
%
% Output
% den is a gird_size(1) *grid_size(2) array contains number density
%%%

global debye_length grid_size number_represent

den = zeros(grid_size(1),grid_size(2));

for i = 1 : n0
    [x,y,w] = weighting(p0(i,1),p0(i,2));
    for j = 1:length(x)
        den(x(j),y(j)) = den(x(j),y(j)) + number_represent * w(j) /debye_length^3;
    end
end
end

% Gauss`s Law correction of E
function E = E_correction(E_0,E_last,den_ref)
%%%
% This function is about using Gauss`s law to correct E, based on
% A. Bruce Langdon,
% On enforcing Gauss' law in electromagnetic particle-in-cell codes,
% Computer Physics Communications,
% Volume 70, Issue 3,
% 1992,
% Pages 447-450,
% ISSN 0010-4655,
% https://doi.org/10.1016/0010-4655(92)90105-8.
% (https://www.sciencedirect.com/science/article/pii/0010465592901058)
%
% place this function after function of solve E and update particle
% position
% Input
% E0 is a structure contains E0.x and E0.y which are gird_size(1) *
% gird_size(2) 2D arrays contains x components and y components of the E
% field
% E_last is a structure contains E0.x and E0.y which are gird_size(1) *
% gird_size(2) 2D arrays contains x components and y components of the E
% field in last step
% den_ref is a gird_size(1) * gird_size(2) 2D array contains charge density
%
% Output
% E is a structure contains E0.x and E0.y which are gird_size(1) *
% gird_size(2) 2D arrays contains x components and y components of the
% corrected E field
%%%

global d debye_length t_unit
[dExdx,~] = gradient(E_last.x,debye_length);
[~,dEydy] = gradient(E_last.y,debye_length);
[dx,dy] = gradient(((dExdx+dEydy)-den_ref),debye_length);
E.x = E_0.x+t_unit*d*dx;
E.y = E_0.y+t_unit*d*dy;
end
E.x = zeros(grid_size(1),grid_size(2));
E.y = zeros(grid_size(1),grid_size(2));
B = zeros(grid_size(1),grid_size(2));
J.x = zeros(grid_size(1),grid_size(2));
J.y = zeros(grid_size(1),grid_size(2));

% Set external field
global E_ext B_ext
E_ext.x = zeros(grid_size(1),grid_size(2));
E_ext.y = zeros(grid_size(1),grid_size(2));
B_ext = zeros(grid_size(1),grid_size(2));

% initialze particle
random_number = rand(particle_number_e,2);
random_number_v = randn(particle_number_e,2)-0.5;
particle_e = [(0.25+0.5*random_number(:,1))*x_total,(0.25+0.5*random_number(:,2))*y_total];
particle_v_e = zeros(particle_number_e,2);
for i = 1:particle_number_e
    particle_v_e(i,1) = v0_e*random_number_v(i,1)/sqrt(random_number_v(i,1)^2+random_number_v(i,2)^2);
    particle_v_e(i,2) = v_drift+v0_e*random_number_v(i,2)/sqrt(random_number_v(i,1)^2+random_number_v(i,2)^2);
end
random_number = randn(particle_number_i,2);
random_number_v = randn(particle_number_i,2)-0.5;
particle_i = [(0.25+0.5*random_number(:,1))*x_total,(0.25+0.5*random_number(:,2))*y_total];
particle_v_i = zeros(particle_number_i,2);
for i = 1:particle_number_i
    particle_v_i(i,1) = v0_i*random_number_v(i,1)/sqrt(random_number_v(i,1)^2+random_number_v(i,2)^2);
    particle_v_i(i,2) = v_drift+v0_i*random_number_v(i,2)/sqrt(random_number_v(i,1)^2+random_number_v(i,2)^2);
end

[x_mesh,y_mesh] = meshgrid(0:debye_length*1e3:(y_total-debye_length)*1e3,0:debye_length*1e3:(x_total-debye_length)*1e3);

%ke = zeros(total_step,1);
total_energy_timestep = zeros(total_step,1);
charge_t = zeros(total_step,1);
charge_t_real = zeros(total_step,1);
charge_count_t = zeros(total_step,1);
total_charge = 0
total_charge_real = 0
total_count = 0
for t = 1:total_step
    clc
    simulation_presentage = 100 * t / total_step
    physical_time = t * t_unit
    particle_number_e
    particle_number_i
    
    total_charge
    total_charge_real
    total_count

    charge_t(t) = total_charge;
    charge_t_real(t) = total_charge_real;
    charge_count_t(t) = total_count;
    E_last = E;
    E = E_solver(E,B,J);

    E.x(1,:) = 0;
    E.x(grid_size(1),:) = 0;
    E.x(:,1) = 0;
    E.x(:,grid_size(2)) = 0;

    E.y(1,:) = 0;
    E.y(grid_size(1),:) = 0;
    E.y(:,1) = 0;
    E.y(:,grid_size(2)) = 0;

    den_e = density(particle_number_e,particle_e);
    den_i = density(particle_number_i,particle_i);

    [particle_e,particle_v_e] = particle_motion(particle_number_e,particle_e,particle_v_e,e_p,m_p_e,E,B);
    [particle_i,particle_v_i] = particle_motion(particle_number_i,particle_i,particle_v_i,e_i,m_p_i,E,B);

    [particle_number_e,particle_e,particle_v_e] = boundary(particle_number_e,particle_e,particle_v_e);
    [particle_number_i,particle_i,particle_v_i] = boundary(particle_number_i,particle_i,particle_v_i);


    % calculate kinetic energy
    % for i = 1 : particle_number_i
    %     ke(t) = ke(t) + 0.5*m_p_i*sqrt(particle_v_i(i,1)^2+particle_v_i(i,2)^2);
    % end
    % for i = 1 : particle_number_e
    %     ke(t) = ke(t) + 0.5*m_p_e*sqrt(particle_v_e(i,1)^2+particle_v_e(i,2)^2);
    % end
    %ke(t)

    % Correct E
    E = E_correction(E,E_last,(den_e * e_p/number_represent + den_i * e_i/number_represent));
    % calculate J
    J_e = J_solver(particle_number_e,particle_e,particle_v_e,e_p);
    J_i = J_solver(particle_number_i,particle_i,particle_v_i,e_i);

    J.x = J_e.x + J_i.x;
    J.y = J_e.y + J_i.y;
    %disp(['Total J_e.x: ', num2str(sum(sum(J_e.x)))]);
    %disp(['Total J_i.x: ', num2str(sum(sum(J_i.x)))]);
    %disp(['Total J_e.y: ', num2str(sum(sum(J_e.y)))]);
    %disp(['Total J_i.y: ', num2str(sum(sum(J_i.y)))]);

    B = B_solver(B,E);


    %inject new particle electron
    inject_n = min(max_injection,particle_number_setted-particle_number_e);
    random_number = rand(inject_n,2);
    particle_e(particle_number_e+1:particle_number_e+inject_n,:) =...
        [(1*random_number(:,1))*x_total,(1*random_number(:,2))*y_total];
    random_number_v = randn(inject_n,2)-0.5;
    for i = 1 : inject_n
        leng = sqrt(random_number_v(i,1)^2+random_number_v(i,2)^2);
        particle_v_e(i+particle_number_e,1) = v0_e*random_number_v(i,1)/leng;
        particle_v_e(i+particle_number_e,2) = v_drift+v0_e*random_number_v(i,2)/leng;
    end
    particle_number_e = particle_number_e + inject_n;

    %inject new particle ion
    inject_n = min(max_injection,particle_number_setted-particle_number_i);
    random_number = rand(inject_n,2);
    particle_i(particle_number_i+1:particle_number_i+inject_n,:) =...
        [(1*random_number(:,1))*x_total,(1*random_number(:,2))*y_total];
    random_number_v = randn(inject_n,2)-0.5;
    for i = 1 : inject_n
        leng = sqrt(random_number_v(i,1)^2+random_number_v(i,2)^2);
        particle_v_i(i+particle_number_i,1) = v0_i*random_number_v(i,1)/leng;
        particle_v_i(i+particle_number_i,2) = v_drift+v0_i*random_number_v(i,2)/leng;
    end
    particle_number_i = particle_number_i + inject_n;


    den_e = density(particle_number_e,particle_e);
    den_i = density(particle_number_i,particle_i);
    cmax = max(max(max(den_i)),max(max(den_e)));
    cmin = min(min(min(den_i)),min(min(den_e)));

    charge_ref = den_ref(E);
    %sum(sum(den_e))*debye_length^3
    %sum(sum(den_i))*debye_length^3
    charge_den = (den_e * e_p + den_i * e_i)/number_represent;

    total_charge = sum(sum(charge_ref))*debye_length^3;
    total_charge_real = sum(sum(charge_den))*debye_length^3;
    total_count = particle_number_i*e_i + particle_number_e*e_p;

    % calculate kinetic energy
    ke = 0;
     for i = 1 : particle_number_i
         ke =ke + 0.5*m_p_i*(particle_v_i(i,1)^2+particle_v_i(i,2)^2);
     end
     for i = 1 : particle_number_e
         ke = ke + 0.5*m_p_e*(particle_v_e(i,1)^2+particle_v_e(i,2)^2);
    end

    total_energy = potential(particle_number_i,particle_number_i,particle_i,particle_i,e_i,e_i)+... 
        potential(particle_number_i,particle_number_e,particle_i,particle_e,e_i,e_p)+... 
        potential(particle_number_e,particle_number_e,particle_e,particle_e,e_p,e_p)+... 
        ke
    total_energy_timestep(t) = total_energy;
end
f = figure(1);
f.Position=[100,100,300,350];
tiledlayout(2,1,"TileSpacing","compact")
nexttile
pcolor(x_mesh,y_mesh,den_e,edgecolor='none');
set(gca,"clim",[cmin,cmax])
title("Electron density")
xlabel("x/mm")
ylabel("y/mm")
nexttile
pcolor(x_mesh,y_mesh,den_i,edgecolor='none');
set(gca,"clim",[cmin,cmax])
title("Ion density")
xlabel("x/mm")
ylabel("y/mm")
cb = colorbar();
clim([cmin,cmax])
cb.Layout.Tile="east";
cb.Label.String = "Number density / m^{-2}";
sgtitle(["Particle distribution at",['t = ',num2str(physical_time),' s']])

cmax = max(max(max(charge_ref)),max(max(charge_den)));
cmin = min(min(min(charge_ref)),min(min(charge_den)));
f2=figure(2);
f2.Position=[400,100,300,350];
tiledlayout(2,1,"TileSpacing","compact");
nexttile
contourf(x_mesh,y_mesh,charge_ref,5);
set(gca,"clim",[cmin,cmax])
title("Deduced charge density")
xlabel("x/mm")
ylabel("y/mm")
nexttile
contourf(x_mesh,y_mesh,charge_den,5);
set(gca,"clim",[cmin,cmax])
title("Charge density")
xlabel("x/mm")
ylabel("y/mm")
cb = colorbar();
clim([cmin,cmax])
cb.Layout.Tile="east";
cb.Label.String = "Charge density / Cm^{-2}";
sgtitle(["Charge distribution at",['t = ',num2str(physical_time),' s']])

figure(3)
plot(1:total_step,charge_t,1:total_step,charge_t_real,1:total_step,charge_count_t)
xlabel("Step")
ylabel("Total charge / C")
title("Total Charge with Step")
legend("Deduced charge","Real Charge")
xlim([1,t+1])

figure(4)
plot(1:total_step,total_energy_timestep)
xlabel("Step")
ylabel("Total energy / J")
title("Total energy with Step")