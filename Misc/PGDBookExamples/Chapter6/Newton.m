function all_FF = Newton(x,t,NL_iter,Max_terms,Max_fp_iter,epsilon,epsilon_tilde,k)
    %Picard iteration for the solution of problem 6.1 with the PGD
    %x,t : space and time discretization
    %NL_iter :  number of Picard iterations
    %Max_terms : maximum number of PGD terms in each increment of the solution
    %Max_fp_iter: maximum number of fixed point iteration in the PGD
    %epsilon: termination criterion for the fixed point
    %epsilon_tilde: termination criterion for the enrichment process
    %k: thermal conductivity
    %
    %Copyright (c) 2013, Francisco Chinesta (Ecole Centrale de Nantes), Roland Keunings (Universite catholique de Louvain), Adrien Leygue(CNRS)
    %Author: Adrien Leygue.
    %All rights reserved.
    %See License file for more details
    
    
    %making sure that x and t are column vectors
    x = x(:);
    t = t(:);
    
    %mesh sizes
    Nx= numel(x);
    Nt = numel(t);
    hx = (x(2)-x(1));
    dt = t(2)-t(1);
    
    %Operators along x
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
    
    %Operators along t
    %mass matrix
    wt = dt*[0.5; ones(Nt-2,1);0.5];
    Mt = spdiags(wt,0,Nt,Nt);
    %time derivative matrix
    D1t = Mt*spdiags([-ones(Nt,1) ones(Nt,1)]/dt,[-1 0],Nt,Nt);
    
    %Cell array that will contain the enrichments computed at each non linear iteration
    all_FF=cell(1,NL_iter);
    
    %Problem definition in matrix operator form
    %Operator
    AA = cell(2,2);
    %RHS
    BB = cell(2,1);
    %A priori known terms,... empty in this case
    GG = cell(2,1);
    %Dirichlet BC
    Dirichlet = cell(2,1);
    %Present solution
    FF = cell(2,1);
    %Mass matrices along x and t
    N_NT = cell(2,1);
    
    
    %initial operator
    AA{1,1} = Mx;
    AA{1,2} = k*Kx;
    AA{2,1} = D1t;
    AA{2,2} = Mt;
    %initial RHS
    BB{1} = Mx*ones(Nx,1);
    BB{2} = Mt*ones(Nt,1);
    
    N_NT{1} = Mx;
    N_NT{2} = Mt;
    
    Dirichlet{1} = [1 Nx];
    Dirichlet{2} = 1;
    
    for newt = 1:NL_iter
        %computation of the increment
        INCREMENT = PGD(AA,BB,GG,N_NT,Dirichlet,Max_terms,Max_fp_iter,epsilon,epsilon_tilde);
        Ninc = size(INCREMENT{1},2);
        %update of the operator from the increment
        if(newt<NL_iter)
            %Update to the linear operator
            AA_INC = cell(2,Ninc);
            for i=(1:Ninc)
                AA_INC{1,i} = 2*bsxfun(@times,Mx,INCREMENT{1}(:,i));
                AA_INC{2,i} = Mt*spdiags(INCREMENT{2}(:,i),0,Nt,Nt);
            end
            AA = [AA AA_INC];
            
            %Update RHS
            %Linear contribution of the increment
            BB{1} = [BB{1} -Mx*INCREMENT{1}  -k*Kx*INCREMENT{1}];
            BB{2} = [BB{2}  D1t*INCREMENT{2}    Mt*INCREMENT{2}];
            %Non-linear contribution
            for i=1:Ninc
                BB{1} = [BB{1} -Mx*bsxfun(@times,[2*FF{1} INCREMENT{1}],INCREMENT{1}(:,i))];
                BB{2} = [BB{2} Mt*bsxfun(@times, [  FF{2} INCREMENT{2}],INCREMENT{2}(:,i))];
            end
        end
        %update solution
        FF = cellfun(@horzcat,FF,INCREMENT,'uniformoutput',false);
        all_FF{newt} = INCREMENT;
    end