%This script calls the function PGD_TransientHeat_TX which implements the
%solution of the transient heat problem 4.1 in order to generate all the
%figures of section 4.4
%
%Copyright (c) 2013, Francisco Chinesta (Ecole Centrale de Nantes), Roland Keunings (Universite catholique de Louvain), Adrien Leygue(CNRS)
%Author: Adrien Leygue.
%All rights reserved.
%See License file for more details


close all
clear all

%domain size & meshes
Lx = 1;
Lt = 0.1;
Nx = 61;
Nt = 151;
x = linspace(0,Lx,Nx)';
t = linspace(0,Lt,Nt)';

%PGD parameters
%Fixed point tolerance
epsilon = 1e-8;
%PGD enrichment tolerance set to 0 to enforce the number of enrichments
epsilon_tilde = 0;
%Maximum dumber of enrichments
Max_terms = 10;
%Maximum number of iterations in the fixed point loop
Max_fp_iter = 50;


%Source term in separated form
FX = ones(size(x));
FT = ones(size(t));

%thermal diffusivity
k=1;
%Reference solution Eq. (4.35)
U_ex = reference_solution_4_35(t,x,k,200);

%PGD solution
[X,T] = PGD_TransientHeat_TX(t,x,Max_terms,Max_fp_iter,epsilon,epsilon_tilde,FT,FX);

%Normalization of the computed enrichment terms
Xnorm = bsxfun(@rdivide,X,sqrt(trapz(x,X.^2)));
Tnorm = bsxfun(@rdivide,T,sqrt(trapz(t,T.^2)));

%Plot X_i
%Figure 4.1
figure;
set(gca,'fontsize',14);
handles=plot(x,Xnorm(:,1:4));
legend(handles,cellfun(@(in) ['X_' num2str(in) '(x)'],num2cell(1:4),'uniformoutput',false))
xlabel('$x$','interpreter','latex')
ylabel('$\frac{X_i(x)}{\|X_i(x)\| \,\,\, \,\,}$','interpreter','latex','fontsize',18)

%Plot T_i
%Figure 4.2
figure;
set(gca,'fontsize',14);
handles=plot(t,Tnorm(:,1:4));
legend(handles,cellfun(@(in) ['T_' num2str(in) '(t)'],num2cell(1:4),'uniformoutput',false))
xlabel('$t$','interpreter','latex')
ylabel('$\frac{T_i(t)}{\|T_i(t)\| \,\,\, \,\,}$','interpreter','latex','fontsize',18)

%Computation of the error berween the reference and PGD solution
E_N = zeros(1,Max_terms);
for j=1:(Max_terms)
    %reconstruction of the PGD solution
    U_pgd = X(:,1:j)*T(:,1:j)';
    E_N(j) = trapz(t,trapz(x,(U_ex - U_pgd).^2));
end

%Figure 4.3
figure;
set(gca,'fontsize',14);
semilogy(1:(Max_terms),E_N);
xlabel('$N$','interpreter','latex')
ylabel('$E(u_{N})$','interpreter','latex')