%This script calls the function PGD_TransientHeat_TXK which implements the
%parametric solution of the transient heat problem 4.1 in order to generate all the
%figures of section 5.1.1
%
%Copyright (c) 2013, Francisco Chinesta (Ecole Centrale de Nantes), Roland Keunings (Universite catholique de Louvain), Adrien Leygue(CNRS)
%Author: Adrien Leygue.
%All rights reserved.
%See License file for more details


close all
clear all

%domain sizes & meshes
Lx = 1;
Lt = 0.1;
Nx = 61;
Nt = 151;
Nk = 101;
x = linspace(0,Lx,Nx)';
t = linspace(0,Lt,Nt)';
k = linspace(1,5,Nk)';

%PGD parameters
%Fixed point tolerance
epsilon = 1e-8;
%PGD enrichment tolerance
epsilon_tilde = 0;
%Maximum dumber of enrichments
Max_terms = 25;
%Maximum number of iterations in the fixed point loop
Max_fp_iter = 50;

%PGD solution
[T,X,K] = PGD_TransientHeat_TXK(t,x,k,Max_terms,Max_fp_iter,epsilon,epsilon_tilde);

%Visualization of the enrichment terms
Xnorm = bsxfun(@rdivide,X,sqrt(trapz(x,X.^2)));
Tnorm = bsxfun(@rdivide,T,sqrt(trapz(t,T.^2)));
Knorm = bsxfun(@rdivide,K,sqrt(trapz(k,K.^2)));

%Plot X_i
%Figure 5.1
figure;
set(gca,'fontsize',14);
handles=plot(x,Xnorm(:,1:4));
legend(handles,cellfun(@(in) ['X_' num2str(in) '(x)'],num2cell(1:4),'uniformoutput',false))
xlabel('$x$','interpreter','latex')
ylabel('$\frac{X_i(x)}{\|X_i(x)\| \,\,\, \,\,}$','interpreter','latex','fontsize',18)

%Plot T_i
%Figure 5.2
figure;
set(gca,'fontsize',14);
handles=plot(t,Tnorm(:,1:4));
legend(handles,cellfun(@(in) ['T_' num2str(in) '(t)'],num2cell(1:4),'uniformoutput',false))
xlabel('$t$','interpreter','latex')
ylabel('$\frac{T_i(t)}{\|T_i(t)\| \,\,\, \,\,}$','interpreter','latex','fontsize',18)

%Plot K_i
%Figure 5.3
figure;
set(gca,'fontsize',14);
handles=plot(k,Knorm(:,1:4));
legend(handles,cellfun(@(in) ['K_' num2str(in) '(k)'],num2cell(1:4),'uniformoutput',false))
xlabel('$t$','interpreter','latex')
ylabel('$\frac{K_i(t)}{\|K_i(t)\| \,\,\, \,\,}$','interpreter','latex','fontsize',18)


%PGD post compression initialization
%cf. example of Chapter 3
AA = cell(3,1);
BB = cell(3,1);
GG = cell(3,1);
N_NT = cell(3,1);
Dirichlet = cell(3,1);
AA{1} = speye(Nt);
AA{2} = speye(Nx);
AA{3} = speye(Nk);
N_NT{1} = speye(Nt);
N_NT{2} = speye(Nx);
N_NT{3} = speye(Nk);

%The RHS is the non-optimal separated representation to compress
BB{1} = T;
BB{2} = X;
BB{3} = K;
%PGD post compression
[FF] = PGD(AA,BB,GG,N_NT,Dirichlet,Max_terms,Max_fp_iter,epsilon,epsilon_tilde);
%New separated representation
T2 = FF{1};
X2 = FF{2};
K2 = FF{3};

%Computation of the error berween the reference and PGD solution for all k
E_N_k = zeros(Max_terms,Nk);
E_N_k2 = zeros(Max_terms,Nk);
%loop over all k
for ik = 1:Nk
    U_ex = reference_solution_5_24(t,x,k(ik),200);
    
    %particularization of the K functions at k(ik)
    Kk = K(ik,:);
    Kk2 = K2(ik,:);
    
    for iN = 1:Max_terms
        %Reconstruction
        U_pgd = bsxfun(@times,X(:,1:iN),Kk(1:iN))*T(:,1:iN)';
        U_pgd2 = bsxfun(@times,X2(:,1:iN),Kk2(1:iN))*T2(:,1:iN)';
        
        E_N_k(iN,ik) = trapz(t,trapz(x,(U_ex - U_pgd).^2));
        E_N_k2(iN,ik) = trapz(t,trapz(x,(U_ex - U_pgd2).^2));
        
    end
end
%integral along k
E_N = trapz(k,E_N_k,2);
E_N2 = trapz(k,E_N_k2,2);

%figure 5.4
figure;
semilogy(1:(Max_terms),E_N,'-',1:(Max_terms),E_N2,'--');
legend('Original solution','Compressed solution')
xlabel('$N$','interpreter','latex')
ylabel('$E(u_{N})$','interpreter','latex')
