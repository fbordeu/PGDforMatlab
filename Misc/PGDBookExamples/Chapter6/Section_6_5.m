%This script calls the function PGD for the problem 6.92 in order to generate
%Figures 6.3 and 6.4
%
%Copyright (c) 2013, Francisco Chinesta (Ecole Centrale de Nantes), Roland Keunings (Universite catholique de Louvain), Adrien Leygue(CNRS)
%Author: Adrien Leygue.
%All rights reserved.
%See License file for more details



close all
clear all
%Steady parametric advection diffusion problem
%with & without upwinding
%Meshes definition
Lx = 2;
Ly = 2;
%Number of dofs in each space coordinate
Nx = 41;
Ny = 41;
%Number of dofs in the velocity coordinates
Nv = 21;
%Number of dofs in the diffusivity coordinate
Nk = 101;
%Meshes
x = linspace(0,Lx,Nx)';
y = linspace(0,Ly,Ny)';
vx = linspace(0.5,1,Nv)';
vy = linspace(0.5,1,Nv)';
%logarithmically spaced grid along k
k = logspace(-3,-1,Nk)';

%x & y mesh grid size
hx = x(2)-x(1);
hy = y(2)-y(1);
%One dimensional operators along x and y
%Identity
Mx = speye(Nx);
My = speye(Ny);

tmp = ones(Nx,1);
%second derivative along x centered scheme
D2x =  spdiags([tmp -2*tmp tmp],[-1 0 1],Nx,Nx) / (hx^2);
D2x(end,end) = D2x(end,end)/2;
%first derivative along x centered scheme
D1x =  spdiags([-tmp  tmp],[-1  1],Nx,Nx) / (2*hx);
%first derivative along y upwind scheme
D1x_up =  spdiags([-tmp  tmp],[-1  0],Nx,Nx) / (hx);

tmp = ones(Ny,1);
%second derivative along y centered scheme
D2y =  spdiags([tmp -2*tmp tmp],[-1 0 1],Ny,Ny) / (hy^2);
D2y(end,end) = D2y(end,end)/2;
%first derivative along y centered scheme
D1y =  spdiags([-tmp  tmp],[-1  1],Ny,Ny) / (2*hy);
%first derivative along y upwind scheme
D1y_up =  spdiags([-tmp  tmp],[-1  0],Ny,Ny) / (hy);

%Identity operators along v and k
Mk = speye(Nk);
Mv = speye(Nv);

%Source term in separated form
fx = 10*exp(- 100*(x-Lx/2).^2);
fy = exp(- 100*(y-Ly/2).^2);
fvx = ones(Nv,1);
fvy = ones(Nv,1);
fk = ones(Nk,1);

%Dims: x, y, vx,vy,k 

%Problem definition in matrix operator form
%Operator
AA = cell(5,4);
%RHS
BB = cell(5,1);
%A priori known terms, to enforce the Dirichlet BC
GG = cell(5,1);
GG{1} = linspace(0,1,Nx)';
GG{2} = linspace(0,1,Ny)';
GG{3} = ones(Nv,1);
GG{4} = ones(Nv,1);
GG{5} = ones(Nk,1);

%Dirichlet BC
Dirichlet = cell(5,1);
%Mass matrice
N_NT = cell(5,1);

%x-diffusion
AA{1,1} = -D2x;
AA{2,1} = My;
AA{3,1} = Mv;
AA{4,1} = Mv;
AA{5,1} = spdiags(k,0,Nk,Nk);
%y-diffusion
AA{1,2} = Mx;
AA{2,2} = -D2y;
AA{3,2} = Mv;
AA{4,2} = Mv;
AA{5,2} = spdiags(k,0,Nk,Nk);
%x-advection, un-stabilized scheme
AA{1,3} = D1x;
AA{2,3} = My;
AA{3,3} = spdiags(vx,0,Nv,Nv);
AA{4,3} = Mv;
AA{5,3} = Mk;
%y-advection, un-stabilized scheme
AA{1,4} = Mx;
AA{2,4} = D1y;
AA{3,4} = Mv;
AA{4,4} = spdiags(vy,0,Nv,Nv);
AA{5,4} = Mk;
%RHS
BB{1} = [fx 0*ones(Nx,1)];
BB{2} = [fy ones(Ny,1)];
BB{3} = [fvx ones(Nv,1)];
BB{4} = [fvy ones(Nv,1)];
BB{5} = [fk ones(Nk,1)];


%Imposition of a Dirichlet BC on the whole boundary
Dirichlet{1} = [1 Nx];
Dirichlet{2} = [1 Ny];

N_NT{1} = Mx;
N_NT{2} = My;
N_NT{3} = Mv;
N_NT{4} = Mv;
N_NT{5} = Mk;

%PGD parameters
%Maximum number of enrichment steps
N_enrich = 200;
%Maximum number of fixed point iterations
N_fp = 55;
%Fixed point iterations termination criterion
epsilon = 1e-8;
%Enrichment loop termination criterion
epsilon_tilde = 1e-8;

%Solution without stabilization, PGD_SYM implements the residual minimization
%form of the PGD
[FF] = PGD_SYM(AA,BB,GG,N_NT,Dirichlet,N_enrich,N_fp,epsilon,epsilon_tilde);
alph = ones(size(FF{1},2),1);

%Figure 6.3
figure;

%reconstruction of U for different cases
%alpha is a vector giving the weight of each term of the separated solution
%after having par'icularized vx, vy and k to a specific value
alpha = FF{3}(1,:).*FF{4}(1,:).*FF{5}(1,:).*alph';
U1 = FF{2}*diag(alpha)*FF{1}';
alpha = FF{3}(1,:).*FF{4}(1,:).*FF{5}(end,:).*alph';
U2 = FF{2}*diag(alpha)*FF{1}';
alpha = FF{3}(1,:).*FF{4}(end,:).*FF{5}(1,:).*alph';
U3 = FF{2}*diag(alpha)*FF{1}';
alpha = FF{3}(end,:).*FF{4}(1,:).*FF{5}(1,:).*alph';
U4 = FF{2}*diag(alpha)*FF{1}';

%values of the isovalues contours
%first solution
%values of the isovalues contours
values = linspace(0,max(U1(:)),10);
subplot(2,2,1),contourf(x,y,U1,values(2:end))
axis equal
set(gca,'fontsize',16);
axis([0 2 0 2]);
xlabel('x','fontsize',16);
ylabel('y','fontsize',16);
set(gca,'xtick',[0 1 2]);
set(gca,'ytick',[0 1 2]);
text(1,0.1,[ '$$v_x=0.5$$, $$v_y=0.5$$' char(10) '$$k=10^{-3}$$'],'interpreter','latex','horizontalalignment','center','verticalalignment','bottom','fontsize',14);

%second solution
values = linspace(0,max(U2(:)),10);
subplot(2,2,2),contourf(x,y,U2,values(2:end))
axis equal
set(gca,'fontsize',16);
axis([0 2 0 2]);
xlabel('x','fontsize',16);
ylabel('y','fontsize',16);
set(gca,'xtick',[0 1 2]);
set(gca,'ytick',[0 1 2]);
text(1,0.1,[ '$$v_x=0.5$$, $$v_y=0.5$$' char(10) '$$k=10^{-1}$$'],'interpreter','latex','horizontalalignment','center','verticalalignment','bottom','fontsize',14);

%third solution
values = linspace(0,max(U3(:)),10);
subplot(2,2,3),contourf(x,y,U3,values(2:end))
axis equal
set(gca,'fontsize',16);
axis([0 2 0 2]);
xlabel('x','fontsize',16);
ylabel('y','fontsize',16);
set(gca,'xtick',[0 1 2]);
set(gca,'ytick',[0 1 2]);
text(1,0.1,[ '$$v_x=0.5$$, $$v_y=1$$' char(10) '$$k=10^{-3}$$'],'interpreter','latex','horizontalalignment','center','verticalalignment','bottom','fontsize',14);


values = linspace(0,max(U4(:)),10);
subplot(2,2,4),contourf(x,y,U4,values(2:end))
axis equal
set(gca,'fontsize',16);
axis([0 2 0 2]);
xlabel('x','fontsize',16);
ylabel('y','fontsize',16);
set(gca,'xtick',[0 1 2]);
set(gca,'ytick',[0 1 2]);
text(1,0.1,[ '$$v_x=1$$, $$v_y=0.5$$' char(10) '$$k=10^{-3}$$'],'interpreter','latex','horizontalalignment','center','verticalalignment','bottom','fontsize',14);



%Modification of the operator to use a stabilized scheme
AA{1,3} = D1x_up;
AA{2,4} = D1y_up;


%Solution with stabilization, PGD_SYM implements the residual minimization
%form of the PGD
[FF] = PGD_SYM(AA,BB,GG,N_NT,Dirichlet,N_enrich,N_fp,epsilon,epsilon_tilde);
alph = ones(size(FF{1},2),1);

%Figure 6.Y
figure;

%reconstruction of U for different cases
%alpha is a vector giving the weight of each term of the separated solution
%after having par'icularized vx, vy and k to a specific value
alpha = FF{3}(1,:).*FF{4}(1,:).*FF{5}(1,:).*alph';
U1 = FF{2}*diag(alpha)*FF{1}';
alpha = FF{3}(1,:).*FF{4}(1,:).*FF{5}(end,:).*alph';
U2 = FF{2}*diag(alpha)*FF{1}';
alpha = FF{3}(1,:).*FF{4}(end,:).*FF{5}(1,:).*alph';
U3 = FF{2}*diag(alpha)*FF{1}';
alpha = FF{3}(end,:).*FF{4}(1,:).*FF{5}(1,:).*alph';
U4 = FF{2}*diag(alpha)*FF{1}';


%first solution
%values of the isovalues contours
values = linspace(0,max(U1(:)),10);
subplot(2,2,1),contourf(x,y,U1,values(2:end))
axis equal
set(gca,'fontsize',16);
axis([0 2 0 2]);
xlabel('x','fontsize',16);
ylabel('y','fontsize',16);
set(gca,'xtick',[0 1 2]);
set(gca,'ytick',[0 1 2]);
text(1,0.1,[ '$$v_x=0.5$$, $$v_y=0.5$$' char(10) '$$k=10^{-3}$$'],'interpreter','latex','horizontalalignment','center','verticalalignment','bottom','fontsize',14);

%second solution
values = linspace(0,max(U2(:)),10);
subplot(2,2,2),contourf(x,y,U2,values(2:end))
axis equal
set(gca,'fontsize',16);
axis([0 2 0 2]);
xlabel('x','fontsize',16);
ylabel('y','fontsize',16);
set(gca,'xtick',[0 1 2]);
set(gca,'ytick',[0 1 2]);
text(1,0.1,[ '$$v_x=0.5$$, $$v_y=0.5$$' char(10) '$$k=10^{-1}$$'],'interpreter','latex','horizontalalignment','center','verticalalignment','bottom','fontsize',14);

%third solution
values = linspace(0,max(U3(:)),10);
subplot(2,2,3),contourf(x,y,U3,values(2:end))
axis equal
set(gca,'fontsize',16);
axis([0 2 0 2]);
xlabel('x','fontsize',16);
ylabel('y','fontsize',16);
set(gca,'xtick',[0 1 2]);
set(gca,'ytick',[0 1 2]);
text(1,0.1,[ '$$v_x=0.5$$, $$v_y=1$$' char(10) '$$k=10^{-3}$$'],'interpreter','latex','horizontalalignment','center','verticalalignment','bottom','fontsize',14);


values = linspace(0,max(U4(:)),10);
subplot(2,2,4),contourf(x,y,U4,values(2:end))
axis equal
set(gca,'fontsize',16);
axis([0 2 0 2]);
xlabel('x','fontsize',16);
ylabel('y','fontsize',16);
set(gca,'xtick',[0 1 2]);
set(gca,'ytick',[0 1 2]);
text(1,0.1,[ '$$v_x=1$$, $$v_y=0.5$$' char(10) '$$k=10^{-3}$$'],'interpreter','latex','horizontalalignment','center','verticalalignment','bottom','fontsize',14);
