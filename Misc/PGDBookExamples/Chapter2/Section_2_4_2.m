%This script calls the function PGD which implements the PGD
%solution of a general problem in order to solve the multidimensional
%Poisson equation and generate all the figures of section 2.4.2
%
%Copyright (c) 2013, Francisco Chinesta (Ecole Centrale de Nantes), Roland Keunings (Universite catholique de Louvain), Adrien Leygue(CNRS)
%Author: Adrien Leygue.
%All rights reserved.
%See License file for more details

%Mesh definition, same for all dimensions
Nx = 101;
x = linspace(-1,1,Nx)';
hx = (x(2)-x(1));

%elemental stiffness matrix
Kex = reshape([1.0./hx,-1.0./hx,-1.0./hx,1.0./hx],[2,2]);
%elemental mass matrix
Mex = reshape([hx.*(1.0./3.0),hx.*(1.0./6.0),hx.*(1.0./6.0),hx.*(1.0./3.0)],[2,2]);

%Mass matrix along each dimension
Mx = sparse(Nx,Nx);
%Stiffness matrix along each dimension
Kx = sparse(Nx,Nx);

%FE assembly of the FE matrices
for el = 1:Nx-1
    Mx([el el+1],[el el+1]) = Mx([el el+1],[el el+1]) + Mex;
    Kx([el el+1],[el el+1]) = Kx([el el+1],[el el+1]) + Kex;
end

%PGD parameters
%Maximum number of enrichment steps
N_enrich = 5;
%Maximum number of fixed point iterations
N_fp = 50;
%Fixed point iterations termination criterion
epsilon = 1e-8;
%Enrichment loop termination criterion
epsilon_tilde = 1e-8;

%Dimensions of the problems: 2D, 5D,...
DD = [2 5 10 20];
%Matrix to store all the computed errors
ERR = zeros(N_enrich+1,numel(DD));

%Loop over the different dimensions
for id = 1:numel(DD)
    D = DD(id);
    
    %Exact solution container
    FF_ex = cell(D,1);
    
    %Container for the differential operator
    AA = cell(D,D);
    %Container for the RHS
    BB = cell(D,1);
    %Container for the a priori known terms of the PGD solution
    GG = cell(D,1);
    %Container to store the nodes subject to a Dirichlet condition
    Dirichlet = cell(D,1);
    %Container to store the mass matrix of each dimension
    N_NT = cell(D,1);
    
    %Construction of the problem and of the manufactured solution along each dimension
    for d = 1:D
        %operator
        for op = 1:D
            if op==d
                AA{d,op} = Kx;
            else
                AA{d,op} = Mx;
            end
        end
        %Dirichlet BC
        Dirichlet{d} = [1 Nx];
        %Exact solution
        FF_ex{d} = [sin(d*pi*x).*x sin((D+1-d)*pi*x).*x.^2];
        %Mass matrix
        N_NT{d} = Mx;
    end
    %Number of terms in the exact solution
    N_ex = size(FF_ex{1},2);
    
    %compute RHS by applying the operators in AA on FF_ex
    for d = 1:D
        BB{d} = zeros(Nx,D*N_ex);
        for op = 1:D
            BB{d}(:,(1:N_ex)+N_ex*(op-1)) = AA{d,op}*FF_ex{d};
        end
    end
    
    %PGD solution
    [FF] = PGD(AA,BB,GG,N_NT,Dirichlet,N_enrich,N_fp,epsilon,epsilon_tilde);
    
    %Compute the quadratic error as a function of the number of enrichment terms
    for i = 0:size(FF{1},2)
        %extract the first i terms of the solution
        SOL = cellfun(@(a) a(:,1:i),FF,'uniformoutput',false);
        %Compute the error
        ERR(i+1,id) = Compute_quadratic_error(FF_ex,SOL,N_NT);
    end
    
end
%Generate figure 2.12
figure;
set(gca,'fontsize',14);
handles=semilogy(0:N_enrich,ERR);
legend(handles,cellfun(@(in) ['$\mathcal{D} = $' num2str(in)],num2cell(DD),'uniformoutput',false),'interpreter','latex');
xlabel('$N$','interpreter','latex')
ylabel('${E_M^\mathcal{D}(u_{M,N}^\mathcal{D}\,)}$','interpreter','latex')


