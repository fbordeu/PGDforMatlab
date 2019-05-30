%This script calls the function PGD_Poisson_2D which implements the
%solution of the Poisson Equation in 2D in order to generate all the
%figures of section 2.1.3
%
%Copyright (c) 2013, Francisco Chinesta (Ecole Centrale de Nantes), Roland Keunings (Universite catholique de Louvain), Adrien Leygue(CNRS)
%Author: Adrien Leygue.
%All rights reserved.
%See License file for more details


close all
clear all
%domain sizes
Lx = 2;
Ly = 1;

%PGD parameters
epsilon = 1e-8;
%PGD enrichment tolerance
epsilon_tilde = 1e-5;
%Maximum number of enrichments
Max_terms = 4;
%Maximum number of iterations in the fixed point loop
Max_fp_iter = 20;

%Mesh
Nx = 101;
Ny = 101;
x = linspace(0,Lx,Nx)';
y = linspace(0,Ly,Ny)';

%PGD solution
[X,Y] = PGD_Poisson_2D(x,y,Max_terms,Max_fp_iter,epsilon,epsilon_tilde);

%Normalization of the computed functions
Xnorm = bsxfun(@rdivide,X,sqrt(trapz(x,X.^2)));
Ynorm = bsxfun(@rdivide,Y,sqrt(trapz(y,Y.^2)));

%Plot X_i
%Figure 2.1
figure;
set(gca,'fontsize',14);
handles=plot(x,Xnorm);
legend(handles,cellfun(@(in) ['X_' num2str(in) '(x)'],num2cell(1:Max_terms),'uniformoutput',false))
xlabel('$x$','interpreter','latex')
ylabel('$\frac{X_i(x)}{\|X_i(x)\| \,\,\, \,\,}$','interpreter','latex','fontsize',18)

%Plot Y_i
%Figure 2.2
figure;
set(gca,'fontsize',14);
handles=plot(y,Ynorm);
legend(handles,cellfun(@(in) ['Y_' num2str(in) '(y)'],num2cell(1:Max_terms),'uniformoutput',false))
xlabel('$y$','interpreter','latex')
ylabel('$\frac{Y_i(y)}{\|Y_i(y)\| \,\,\, \,\,}$','interpreter','latex','fontsize',18)

%Plot reconstructed solution
%Figure 2.3
figure;
set(gca,'fontsize',14);
hold on
plot3(x,zeros(size(x)),Xnorm/50,'linewidth',2);
plot3(zeros(size(y)),y,Ynorm/50,'linewidth',2);
surf(x,y,Y*X');
colormap('cool');
shading interp
xlabel('x');
ylabel('y');
zlabel('U_4(x,y)')
view(3);



%Computation of the L2 Error as a function of the discretization & enrichment
%
%values of the tolerances for the convergence
%fixed point iterations tolerance
epsilon = 1e-8;
%PGD enrichment tolerance set to zero to enforce the number of enrichments
epsilon_tilde = 0;
%Maximum dumber of enrichments
Max_terms = 12;
%Maximum number of iterations in the fixed point loop
Max_fp_iter = 20;

%Different values of M (# of grid points per dimension)
M = [11 21 41 61 81 101 121];

%Container for the error E_MN(i,j) contains the error for M(i) grid points and j enrichment terms
E_MN = zeros(numel(M),Max_terms);
%Loop over the discretization
for i=1:numel(M)
    %Mesh
    Nx = M(i);
    Ny = M(i);
    x = linspace(0,Lx,Nx)';
    y = linspace(0,Ly,Ny)';
    
    %PGD solution
    [X,Y] = PGD_Poisson_2D(x,y,Max_terms,Max_fp_iter,epsilon,epsilon_tilde);
    
    %Analytical solution
    U_exact = reference_solution_2_28(x,y,50);
    
    %Error computation
    for j=1:Max_terms
        U_pgd = Y(:,1:j)*X(:,1:j)';
        E_MN(i,j) = trapz(x,trapz(y,(U_exact - U_pgd).^2,1),2);
    end
    
    %Error maps for M = M(3)
    %Figure 2.4
    if(i==3)
        figure;
        set(gca,'fontsize',14);
        
        %N=1
        subplot(2,2,1),
        surf(x,y,U_exact-Y(:,1)*X(:,1)');
        title('N=1','interpreter','latex');
        xlabel('x','interpreter','latex');
        ylabel('y','interpreter','latex');
        zlabel('$u_{M,N}-u_{ex}$','interpreter','latex','fontsize',18);
        colormap('cool');
        
        %N=2
        subplot(2,2,2),
        surf(x,y,U_exact-Y(:,1:2)*X(:,1:2)');
        title('N=2','interpreter','latex');
        xlabel('x','interpreter','latex');
        ylabel('y','interpreter','latex'),
        zlabel('$u_{M,N}-u_{ex}$','interpreter','latex','fontsize',18);
        colormap('cool');
        
        %N=3
        subplot(2,2,3),
        surf(x,y,U_exact-Y(:,1:3)*X(:,1:3)');
        title('N=3','interpreter','latex');
        xlabel('x','interpreter','latex');
        ylabel('y','interpreter','latex');
        zlabel('$u_{M,N}-u_{ex}$','interpreter','latex','fontsize',18);
        colormap('cool');
        
        %N=10
        subplot(2,2,4),
        surf(x,y,U_exact-Y*X');
        title('N=10','interpreter','latex');
        xlabel('x','interpreter','latex');
        ylabel('y','interpreter','latex');
        zlabel('$u_{M,N}-u_{ex}$','interpreter','latex','fontsize',18);
        colormap('cool');
        
        set(gcf,'renderer','painters')
    end
end


%Error as a function of N for different values of M
%Figure 2.5
figure;
set(gca,'fontsize',14);
handles = semilogy(1:Max_terms,E_MN);
legend(handles,cellfun(@(in) ['M=' num2str(in)],num2cell(M),'uniformoutput',false))
set(gca,'xtick',[0 2 4 6 8 10]);
set(gca,'xlim',[0 Max_terms+1]);
xlabel('$N$','interpreter','latex')
ylabel('$E_M(u_{M,N})$','interpreter','latex','fontsize',18)

%Error as a function of M for different numbers of enrichment
%Figure 2.6
figure
set(gca,'fontsize',14);
Ni = [1 2 3 4 5 6 10];
handles=semilogy(M,E_MN(:,Ni));
legend(handles,cellfun(@(in) ['N=' num2str(in)],num2cell(Ni),'uniformoutput',false))
set(gca,'xtick',M);
xlabel('$M$','interpreter','latex')
ylabel('$E_M(u_{M,N})$','interpreter','latex','fontsize',18)



%Figure 2.7
%Error as a function of the discretization & number of modes & mesh refinment
%PGD parameters
epsilon = 1e-8;
%PGD enrichment tolerance
epsilon_tilde = 0;
%Maximum dumber of enrichments for the different meshes
Max_terms_coarse = 10;
Max_terms_fine = 10;

%Maximum number of iterations in the fixed point loop
Max_fp_iter = 20;

%container for the error
E = zeros(1,(Max_terms_coarse+Max_terms_fine));

%Coarse mesh
Nx_coarse = 21;
Ny_coarse = 21;
x_coarse = linspace(0,Lx,Nx_coarse)';
y_coarse = linspace(0,Ly,Ny_coarse)';

%Fine mesh
Nx_fine = 41;
Ny_fine = 41;
x_fine = linspace(0,Lx,Nx_fine)';
y_fine = linspace(0,Ly,Ny_fine)';


%PGD solution on coarse mesh
[X_coarse,Y_coarse] = PGD_Poisson_2D(x_coarse,y_coarse,Max_terms_coarse,Max_fp_iter,epsilon,epsilon_tilde);

%Interpolation of the functions on the coarse mesh on the fine mesh
X_fine_0 = interp1(x_coarse,X_coarse,x_fine,'spline');
Y_fine_0 = interp1(y_coarse,Y_coarse,y_fine,'spline');

%PGD solution on the fine mesh, computation of additional enrichments
[X_fine,Y_fine] = PGD_Poisson_2D(x_fine,y_fine,Max_terms_fine,Max_fp_iter,epsilon,epsilon_tilde,X_fine_0,Y_fine_0);


%Error with the coarse enrichments
U_exact = reference_solution_2_28(x_coarse,y_coarse,50);
for j= (1):(Max_terms_coarse)
    U_pgd = Y_coarse(:,1:j)*X_coarse(:,1:j)';
    E(j) = trapz(x_coarse,trapz(y_coarse,(U_exact - U_pgd).^2,1),2);
end

%Error with the coarse enrichments + fine enrichments
U_exact = reference_solution_2_28(x_fine,y_fine,50);
for j= (Max_terms_coarse+1):(Max_terms_coarse+Max_terms_fine)
    U_pgd = Y_fine(:,1:j)*X_fine(:,1:j)';
    E(j) = trapz(x_fine,trapz(y_fine,(U_exact - U_pgd).^2,1),2);
end


%figure 2.7
figure
set(gca,'fontsize',14);
semilogy(1:(Max_terms_coarse+Max_terms_fine),E);
h = line('xdata',[1 1]*(Max_terms_coarse+0.5),'ydata',get(gca,'ylim'));
set(h,'linestyle',':','color','black');
text(Max_terms_coarse+0.5,2e-7,'$M=21\, \leftarrow \quad$','interpreter','latex','horizontalalignment','right');
text(Max_terms_coarse+0.5,2e-7,'$\quad \rightarrow \, M=41$','interpreter','latex','horizontalalignment','left');
xlabel('$N$','interpreter','latex')
ylabel('$E_M(u_{M,N})$','interpreter','latex','fontsize',18)



