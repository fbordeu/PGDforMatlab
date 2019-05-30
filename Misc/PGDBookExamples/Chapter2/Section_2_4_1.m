%This script calls the function PGD_Poisson_2D_FEM which implements the
%solution of the Poisson Equation in 2D in order to generate all the
%figures of section 2.4.1
%
%Copyright (c) 2013, Francisco Chinesta (Ecole Centrale de Nantes), Roland Keunings (Universite catholique de Louvain), Adrien Leygue(CNRS)
%Author: Adrien Leygue.
%All rights reserved.
%See License file for more details


close all
clear all

%domain dimensions
Lx = 2;
Ly = 1;

%PGD parameters
epsilon = 1e-8;
%PGD enrichment tolerance set to zero to enforce the number of enrichments
epsilon_tilde = 0;
%Maximum number of enrichments
Max_terms = 10;
%Maximum number of iterations in the fixed point loop
Max_fp_iter = 20;

%Mesh
Nx = 41;
Ny = 41;
x = linspace(0,Lx,Nx)';
y = linspace(0,Ly,Ny)';

%A priori terms imposing the Dirichlet BC (Eq. (2.58))
GX = (Lx-x)/Lx;
GY = y.*(Ly-y);
NG= size(GX,2);

%Separated representation of the source term (Eq. (2.59))
FX = -5*exp(-10*(x-Lx/2).^2);
FY = exp(-10*(y-Ly/2).^2);

%Neumann BC on upper boundary (Eq. (2.57), last line)
q = -1;

%2D Finite element reference solution
U_fem = reference_solution_2_56_FEM(x,y,GX,GY,q,FX,FY);

%PGD solution
[X,Y] = PGD_Poisson_2D_FEM(x,y,Max_terms,Max_fp_iter,epsilon,epsilon_tilde,GX,GY,q,FX,FY);
%NB: The first columns of X and Y are GX and GY respectively.

%Normalization of the first four computed enrichments
Xnorm = bsxfun(@rdivide,X(:,1:5),sqrt(trapz(x,X(:,1:5).^2)));
Ynorm = bsxfun(@rdivide,Y(:,1:5),sqrt(trapz(y,Y(:,1:5).^2)));

%Plot G_1^x, X_i
%Figure 2.8
figure;
set(gca,'fontsize',14);
handles=plot(x,Xnorm);
legend(handles,[{'G_1^x(x)'} cellfun(@(in) ['X_' num2str(in) '(x)'],num2cell(1:4),'uniformoutput',false)])
xlabel('$x$','interpreter','latex')
ylabel('$\frac{X_i(x)}{\|X_i(x)\| \,\,\, \,\,}$','interpreter','latex','fontsize',18)

%Plot G_1^y, Y_i
%Figure 2.9
figure;
set(gca,'fontsize',14);
handles=plot(y,Ynorm);
legend(handles,[{'G_1^y(y)'} cellfun(@(in) ['Y_' num2str(in) '(y)'],num2cell(1:4),'uniformoutput',false)])
xlabel('$y$','interpreter','latex')
ylabel('$\frac{Y_i(y)}{\|Y_i(y)\| \,\,\, \,\,}$','interpreter','latex','fontsize',18)

%Plot of difference with FE solution for different number of enrichments
%Figure 2.10
figure;
set(gca,'fontsize',14)
subplot(2,2,1),
surf(x,y,U_fem-Y(:,1)*X(:,1)');
title('N=0','interpreter','latex');
xlabel('x','interpreter','latex');
ylabel('y','interpreter','latex');
zlabel('$u_{\mathrm{FE} , M} - u_{M,N}$','fontsize',18,'interpreter','latex');
colormap('cool');

subplot(2,2,2),
surf(x,y,U_fem-Y(:,1:2)*X(:,1:2)');
title('N=2','interpreter','latex');
xlabel('x','interpreter','latex');
ylabel('y','interpreter','latex');
zlabel('$u_{\mathrm{FE} , M} - u_{M,N}$','fontsize',18,'interpreter','latex');
colormap('cool');

subplot(2,2,3),
surf(x,y,U_fem-Y(:,1:5)*X(:,1:5)');
title('N=4','interpreter','latex');
xlabel('x','interpreter','latex');
ylabel('y','interpreter','latex');
zlabel('$u_{\mathrm{FE} , M} - u_{M,N}$','fontsize',18,'interpreter','latex');
colormap('cool');

subplot(2,2,4),
surf(x,y,U_fem-Y*X');
title('N=10','interpreter','latex');
xlabel('x','interpreter','latex');
ylabel('y','interpreter','latex');
zlabel('$u_{\mathrm{FE} , M} - u_{M,N}$','fontsize',18,'interpreter','latex');
set(gcf,'renderer','painters')
colormap('cool');


%Computation of the  error between PGD and FEM solution
E_N = zeros(1,Max_terms+size(GX,2));
for j=1:(Max_terms+size(GX,2))
    U_pgd = Y(:,1:j)*X(:,1:j)';
    E_N(j) = trapz(x,trapz(y,(U_fem - U_pgd).^2,1),2);
end

%Figure 2.11
figure;
set(gca,'fontsize',14);
semilogy(0:Max_terms,E_N)
xlabel('$N$','interpreter','latex')
ylabel('$E_M(u_{M,N})$','interpreter','latex','fontsize',18)
