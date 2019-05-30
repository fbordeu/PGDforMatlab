%Code to reproduce the results of Section 1.2.3 and produce Figures 1.1 to
%1.4
%
%Copyright (c) 2013, Francisco Chinesta (Ecole Centrale de Nantes), Roland Keunings (Universite catholique de Louvain), Adrien Leygue(CNRS)
%Author: Adrien Leygue.
%All rights reserved.
%See License file for more details

close all
clear all;

%Parameter definition
%thermal diffusivity
lambda = 0.01;
%Number of eigenvectors to keep in the reduced model
N_eig = 3;

%mesh definition in space & time
Nx = 51;
Nt = 31;
x = linspace(0,1,Nx);
t = linspace(0,30,Nt);
h = x(2)-x(1);
dt = t(2)-t(1);
%coarse mesh indices
coarse_x = 6:5:Nx;
coarse_t = [1     6    11    16    21        31];


%Definition of the FD differentiation matrix, 2nd order centered scheme
e = ones(Nx,1);
D2 = spdiags([e -2*e e]/(h^2),[-1 0 1],Nx,Nx);

%special treatment of the first point to account for the flux boundary
%condition
% u(h) is approximated by u(0) + h* u'(0) + (1/2)*h^2* u''(0)
%since u'(0) = q
%u''(0) is approximated by (2/h^2)*(u(h) - u(0)) -2*q/h
D2(1,[1 2]) = 2*[-1 1]/h^2;

%special treatment of the last point to account for the flux boundary
%condition
% u(1-h) is approximated by u(1) - h*u'(1) + (1/2)*h^2* u''(1)
% since u'(1) = 0
% u''(1) is approximated by (2/h^2)*(u(1-h) - u(1))
D2(end,[Nx-1 Nx]) = 2*[1 -1]/h^2;

%The discrete system writes
%d/dt (U)  = lambda*D2*U + lambda*(2/h)*q(t)*B
%where B is a vector of size Nx of zeros excepted on the first line where
%it is 1

%An implicit euler discretization yields
% (U(t+1)-U(t)) / dt = lambda*D2*U(t+1) + lambda*(2/h)*q(t+1)*B;
%(I - dt*lambda*D2)*U(t+1) = I*U(t) + lambda*dt*(2/h)*q(t+1)*B;
%where I is the identity matrix.
I = speye(Nx,Nx);
%to conform with Eq. 1.12 we divide the RHS & LHS by (lambda*dt*2/h)
%K*U(t+1) = M*U(t) + q(t);
%where
K = (I - dt*lambda*D2)/(lambda*dt*2/h);
M = I/(lambda*dt*2/h);
%source term definition through a local function
q =  @(arg_t) [double(arg_t<=10); zeros(Nx-1,1)];

%The columns of U will contain the value of U at different times
%the initial condition is homogeneous
U = zeros(Nx,Nt);

%solution: Nt-1 time steps
for i=2:Nt
    U(:,i) = K\(M*U(:,i-1) + q(t(i)));
end

%Plot of the reference solution Figure 1.1
figure;
hold on
set(gca,'fontsize',14);
handles = plot(x,U(:,coarse_t));
legend(handles,cellfun(@(in) ['t=' num2str(in)],num2cell(coarse_t),'uniformoutput',false))
xlabel('x');
ylabel('Temperature');
title('Reference solution');
axis([0 1 -0.05 0.4]);

%construction of the matrix c
%Q is U(:,2:end) as the first column uf U is zero
c = U(:,2:end)*U(:,2:end)';
%force symmetry. Non-symmetry arises from floating point arithmetics in the
%construction of c
c = 0.5*(c+c');
%computation of the eigenvectors & eigenvalues of c
[PHI,ALPHA] = eig(c);
%selection only the diagonal elements of the eigenvalues matrix
ALPHA = diag(ALPHA)';
%reverse the ordering of the eigen-vectors/-values to have them sorted by
%decreasing eigenvalue
ALPHA = fliplr(ALPHA);
PHI = fliplr(PHI);

%select only the eigenvectors of the dominant eigenvalues
ALPHA = ALPHA(1:N_eig);
B = PHI(:,1:N_eig);

%Figure 1.2
disp(['Dominant eigenvalues: ' num2str(ALPHA(1:N_eig))]);
figure;
hold on;
set(gca,'fontsize',14);
handles=plot(x,B);
legend(handles,cellfun(@(in) ['\Phi_' num2str(in)],num2cell(1:N_eig),'uniformoutput',false))
xlabel('x');
ylabel('\Phi_i');
title('Dominant eigenvectors');

%Construction of the reduced order system
Kred = B'*K*B;
Mred = B'*M*B;

%The columns of zeta will contain the value of reduced model coefficients at different times
%the initial condition is homogeneous
zeta = zeros(N_eig,Nt);

%Solution using the Reduced Order Model: Nt-1 time steps
for i=2:Nt
    zeta(:,i) = Kred\(Mred*zeta(:,i-1) + B'*q(t(i)));
end

%Figure 1.3
%ROM reconstruction & comparison with the reference solution
figure;
hold on
set(gca,'fontsize',14);
handles = plot(x,U(:,coarse_t));
legend(handles,cellfun(@(in) ['t=' num2str(in)],num2cell(coarse_t),'uniformoutput',false))
plot(x(coarse_x),B(coarse_x,:)*zeta(:,coarse_t),'o');
xlabel('x');
ylabel('Temperature');
title('reference solution (line) & ROM predictions (o)');
axis([0 1 -0.05 0.4]);

%NEW source term definition through a local function
q =  @(arg_t) [double(arg_t<=10)*arg_t/10 + double(arg_t>10)*(arg_t-30)/10      ; zeros(Nx-1,1)];

%initialization of the full model
U = zeros(Nx,Nt);
%initialization of the reduced model
zeta = zeros(N_eig,Nt);

%solution of both models: Nt-1 time steps
for i=2:Nt
    U(:,i) = K\(M*U(:,i-1) + q(t(i)));
    zeta(:,i) = Kred\(Mred*zeta(:,i-1) + B'*q(t(i)));
end

%Figure 1.4
%ROM reconstruction & comparison with the reference solution
figure;
hold on
set(gca,'fontsize',14);
handles=plot(x,U(:,coarse_t));
legend(handles,cellfun(@(in) ['t=' num2str(in)],num2cell(coarse_t),'uniformoutput',false))
plot(x(coarse_x),B(coarse_x,:)*zeta(:,coarse_t),'o');
xlabel('x');
ylabel('Temperature');
title('New problem FD solution (line) & ROM predictions (o)');