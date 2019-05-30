function all_FF = Picard(x,t,NL_iter,Max_terms,Max_fp_iter,epsilon,epsilon_tilde,k)
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
    wt = dt*[0.5; ones(Nt-2,1);0.5];
    %mass matrix
    Mt = spdiags(wt,0,Nt,Nt);
    %differentiation matrix
    D1t = Mt*spdiags([-ones(Nt,1) ones(Nt,1)]/dt,[-1 0],Nt,Nt);
    
    %Cell array that will contain the enrichments computed at each non linear iteration
    all_FF=cell(1,NL_iter);
    
    %Problem definition in matrix operator format
    %Operators
    AA = cell(2,2);
    %Initial RHS
    BB0 = cell(2,1);
    %Mass matrices along x and t
    N_NT = cell(2,1);
    %Dirichlet BC
    Dirichlet = cell(2,1);
    %Current solution
    FF = cell(2,1);
    
    Dirichlet{1} = [1 Nx];
    Dirichlet{2} = 1;
    
    BB0{1} = Mx*ones(Nx,1);
    BB0{2} = Mt*ones(Nt,1);
    %RHS of the original problem
    BB = BB0;
    
    AA{1,1} = Mx;
    AA{1,2} = k*Kx;
    AA{2,1} = D1t;
    AA{2,2} = Mt;
    N_NT{1} = Mx;
    N_NT{2} = Mt;
    
    Nterms=0;
    for pic = 1:NL_iter
        %the a prori known modes describe the solution computed up to this
        %point
        GG = FF;
        FF = PGD(AA,BB,GG,N_NT,Dirichlet,Max_terms+Nterms,Max_fp_iter,epsilon,epsilon_tilde);
        %Update of the RHS:
        %original RHS
        BB = BB0;
        %contribution of the just-computed increment through the non-linear
        %term
        Nterms = size(FF{1},2);
        if pic<NL_iter
            %Computation of the quadratic term on the RHS from FF
            for i=1:Nterms
                BB{1} = [BB{1} -Mx*bsxfun(@times,FF{1},FF{1}(:,i))];
                BB{2} = [BB{2} Mt*bsxfun(@times,FF{2},FF{2}(:,i))];
            end
        end
        %we extract the new enrichments computed in the last call to PGD
        N_new_terms = Nterms - size(GG{1},2);
        %we store them for future use
        all_FF{pic} = cellfun(@(in1) in1(:,(Nterms-N_new_terms+1):Nterms),FF,'uniformoutput',false);
    end
end