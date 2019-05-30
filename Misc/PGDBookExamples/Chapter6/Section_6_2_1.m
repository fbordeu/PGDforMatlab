%This script calls the function PGD for the problem 6.28 in order to generate
%Figure 6.2
%
%Copyright (c) 2013, Francisco Chinesta (Ecole Centrale de Nantes), Roland Keunings (Universite catholique de Louvain), Adrien Leygue(CNRS)
%Author: Adrien Leygue.
%All rights reserved.
%See License file for more details


%Steady advection diffusion problem
%with & without upwinding
close all
clear all

%Meshes definition
Lx = 2;
Ly = 2;
Nx = 41;
Ny = 41;
x = linspace(0,Lx,Nx)';
y = linspace(0,Ly,Ny)';
hx = x(2)-x(1);
hy = y(2)-y(1);

%One dimensional operators along x and y
%Identity
Mx = speye(Nx);
My = speye(Ny);

%Second and first derivatives
tmp = ones(Nx,1);
%second derivative along x centered scheme
D2x =  spdiags([tmp -2*tmp tmp],[-1 0 1],Nx,Nx) / (hx^2);
%first derivative along x centered scheme (i.e. without convective stabilization)
D1x =  spdiags([-tmp  tmp],[-1  1],Nx,Nx) / (2*hx);
%first derivative along y upwind scheme (i.e. with convective stabilization)
D1x_up =  spdiags([-tmp  tmp],[-1  0],Nx,Nx) / (hx);

tmp = ones(Ny,1);
%second derivative along y centered scheme
D2y =  spdiags([tmp -2*tmp tmp],[-1 0 1],Ny,Ny) / (hy^2);
%first derivative along y centered scheme
D1y =  spdiags([-tmp  tmp],[-1  1],Ny,Ny) / (2*hy);
%first derivative along y upwind scheme
D1y_up =  spdiags([-tmp  tmp],[-1  0],Ny,Ny) / (hy);

%Physical parameters
%Diffusivity
k = 1e-3;
%Velocity
vx = 1;
vy = 1;
%Source term
fx = 10*exp(- 100*(x-Lx/2).^2);
fy = exp(- 100*(y-Ly/2).^2);

%Problem definition in matrix operator form
%Operator
AA = cell(2,4);
%RHS
BB = cell(2,1);
%A priori known terms, to enforce the Dirichlet BC
GG = cell(2,1);
%Dirichlet BC
Dirichlet = cell(2,1);
%Mass matrices along x and y
N_NT = cell(2,1);

%x-diffusion
AA{1,1} = -k*D2x;
AA{2,1} = My;
%y-diffusion
AA{1,2} = k*Mx;
AA{2,2} = -D2y;
%x-advection, centered scheme
AA{1,3} = vx*D1x;
AA{2,3} = My;
%y-advection, centered scheme
AA{1,4} = vy*Mx;
AA{2,4} = D1y;

%RHS
BB{1} = Mx*fx;
BB{2} = My*fy;

%A priori known term to enforce the Dirichlet BC
GG{1} = linspace(0,1,Nx)';
GG{2} = linspace(0,1,Ny)';

%Imposition of a Dirichlet BC on the whole boundary
Dirichlet{1} = [1 Nx];
Dirichlet{2} = [1 Ny];

N_NT{1} = Mx;
N_NT{2} = My;

%PGD parameters
%Maximum number of enrichment steps
N_enrich = 40;
%Maximum number of fixed point iterations
N_fp = 50;
%Fixed point iterations termination criterion
epsilon = 1e-8;
%Enrichment loop termination criterion
epsilon_tilde = 1e-8;

%Solution without upwinding
[FF] = PGD(AA,BB,GG,N_NT,Dirichlet,N_enrich,N_fp,epsilon,epsilon_tilde);

%Modification of AA to account for convective stabilization
%x-advection, centered scheme
AA{1,3} = vx*D1x_up;
AA{2,3} = My;
%y-advection, centered scheme
AA{1,4} = vy*Mx;
AA{2,4} = D1y_up;
%Solution with convective stabilization
[FF_up] = PGD(AA,BB,GG,N_NT,Dirichlet,N_enrich,N_fp,epsilon,epsilon_tilde);


%figure 6.2
figure;

%left part, without stabilization
U = FF{2}*FF{1}';
subplot(1,2,1),contourf(x,y,U)
axis equal
set(gca,'fontsize',16);
axis([0 2 0 2]);
colormap(flipud(colormap));
xlabel('x','fontsize',16);
ylabel('y','fontsize',16);
title('Without stabilization')
set(gca,'xtick',[0 1 2]);
set(gca,'ytick',[0 1 2]);

%right part, with stabilization
U_up = FF_up{2}*FF_up{1}';
subplot(1,2,2),contourf(x,y,U_up,[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
axis equal
set(gca,'fontsize',16);
axis([0 2 0 2]);
colormap(flipud(colormap));
xlabel('x','fontsize',16);
ylabel('y','fontsize',16);
title('With stabilization')
set(gca,'xtick',[0 1 2]);
set(gca,'ytick',[0 1 2]);
