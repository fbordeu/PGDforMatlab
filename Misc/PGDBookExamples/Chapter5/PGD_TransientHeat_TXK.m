function [T,X,K] = PGD_TransientHeat_TXK(t,x,k,Max_terms,Max_fp_iter,epsilon,epsilon_tilde)
    %PGD solution of the parametric problem 5.23 with homogeneous Dirichlet
    %boundary conditions and homogeneous initial condition
    %For this example, we use simple second order FD with trapezoidal
    %integration rule along the space dimension, an implicit euler scheme
    %along the time dimension and colocation aling the parametric dimension
    %Outputs:
    %       T,X,K :  computed PGD solution. The ith column of X (resp. T,K) contains the nodal values of X_i (resp. T_i,K_i)
    %Inputs:
    %       t,x,k: uniform 1D grids used along each dimension
    %       Max_terms: Maximum number of enrichments
    %       Max_fp_iter: Maximum number of fixed point iterations
    %       epsilon: termination criterion for the fixed point (Eq. (2.8))
    %       epsilon_tilde: termination criterion used for the enrichment process (Eq. (2.26))
    %
    %Copyright (c) 2013, Francisco Chinesta (Ecole Centrale de Nantes), Roland Keunings (Universite catholique de Louvain), Adrien Leygue(CNRS)
    %Author: Adrien Leygue.
    %All rights reserved.
    %See License file for more details
    
    
    %source term magnitude
    f=1;
    
    %Mesh definition for each dimension
    Nt = numel(t);
    Nx = numel(x);
    Nk = numel(k);
    %reshape to ensure that x & y are column vectors
    t = t(:);
    x = x(:);
    k = k(:);
    %Mesh size
    dt = (t(2)-t(1));
    hx = (x(2)-x(1));
    
    %Discrete operators along each dimension
    %Finite differences matrices in time
    %Identity matrix
    It = speye(Nt);
    %time differentiation matrix
    D1t = spdiags([-ones(size(t)) ones(size(t))]/dt,[-1 0],Nt,Nt);
    D1t(1,:) = 0;
    
    %Finite differences matrices in x
    %Identity matrix
    Ix = speye(Nx);
    %Finite Differences second order differentiation matrix
    D2x = spdiags([ones(Nx,1)  -2*ones(Nx,1)  ones(Nx,1)],[-1 0 1],Nx,Nx)/hx^2;
    
    X = zeros(Nx,0);
    T = zeros(Nt,0);
    K = zeros(Nk,0);
    
    %main enrichment loop
    for term=1:Max_terms
        %initialization of the fixed point loop
        St = randn(Nt,1);
        Sx = randn(Nx,1);
        Sk = randn(Nk,1);
        
        %Satisfaction of the homogeneous Dirichlet Boundary conditions for the
        %enrichments
        St(1) = 0;
        Sx(1) = 0;
        Sx(end)=0;
        
        %precomputation of some scalar quantities to i nitialize the fixed
        %point loop
        s1 = trapz(t,St.^2);
        s2 = trapz(t,St.*(D1t*St));
        s3 = trapz(t,St);
        %In case this is not the first enrichment
        if (term>1)
            %RHS coefficients
            s4_i = trapz(t,bsxfun(@times,St,D1t*T));
            s5_i = trapz(t,bsxfun(@times,St,T));
        end
        
        %fixed point iterations
        for iter=1:Max_fp_iter
            %Store the old values of Sx, St & Sk for later comparison
            St_old = St;
            Sx_old = Sx;
            Sk_old = Sk;
            
            %Solve for Sx
            %construction of the boundary value problem along x
            %coefficients
            w1 = trapz(k,Sk.^2);
            w2 = trapz(k,k.*(Sk.^2));
            w3 = trapz(k,Sk);
            %Construction of the RHS
            RHS = w3*s3*f*ones(size(x));
            
            %In case this is not the first enrichment, previous terms are added
            %to the RHS
            if (term>1)
                %RHS coefficients
                w4_i = trapz(k,bsxfun(@times,Sk,K));
                w5_i = trapz(k,bsxfun(@times,Sk.*k,K));
                
                RHS = RHS - X*(s4_i.*w4_i)' + (D2x*X)*(s5_i.*w5_i)';
            end
            
            %construction of the FD boundary value problem
            A = (w1*s2)*Ix -(w2*s1)*D2x;
            %solution with homogeneous boundary conditions at both ends
            Sx(2:end-1) = A(2:end-1,2:end-1)\RHS(2:end-1);
            
            %Solve for St
            %construction of the initial value problem along t
            %coefficients
            r1 = trapz(x,Sx.^2);
            r2 = trapz(x,Sx.*(D2x*Sx));
            r3 = trapz(x,Sx);
            
            %Source term coefficient
            
            %Construction of the RHS
            RHS = w3*r3*f*ones(size(t));
            %In case this is not the first enrichment, previous terms are added
            %to the RHS
            if (term>1)
                %RHS coefficients
                r4_i = trapz(x,bsxfun(@times,Sx,D2x*X));
                r5_i = trapz(x,bsxfun(@times,Sx,X));
                
                RHS = RHS -(D1t*T)*(w4_i.*r5_i)' + T*(w5_i.*r4_i)';
            end
            
            %construction of the initial value problem
            A = (w1*r1)*D1t - (w2*r2)*It;
            %solution with homogeneous initial condition.
            St(2:end) = A(2:end,2:end)\RHS(2:end);
            
            %Solve for Sk
            %construction of the initial value problem along t
            %coefficients
            s1 = trapz(t,St.^2);
            s2 = trapz(t,St.*(D1t*St));
            s3 = trapz(t,St);
            
            %Construction of the RHS
            RHS = s3*r3*f*ones(Nk,1);
            %In case this is not the first enrichment, previous terms are added
            %to the RHS
            if (term>1)
                %RHS coefficients
                s4_i = trapz(t,bsxfun(@times,St,D1t*T));
                s5_i = trapz(t,bsxfun(@times,St,T));
                
                RHS = RHS -K*(s4_i.*r5_i)' +  bsxfun(@times,K,k)*(r4_i.*s5_i)';
            end
            %construction of the algebraic problem
            diagA = (r1*s2-r2*s1*k);
            %solution.
            Sk = RHS./diagA;
            
            %Norm of the difference between the 2 fixed point iterations (E. (2.8))
            S_difference = sqrt(trapz(x,Sx.^2)*trapz(t,St.*St)*trapz(k,Sk.^2) + trapz(x,Sx_old.^2)*trapz(t,St_old.^2)*trapz(k,Sk_old.^2) - 2*trapz(x,Sx.*Sx_old)*trapz(t,St.*St_old)*trapz(k,Sk.*Sk_old));
            %fixed point exit test
            if(S_difference < epsilon), break; end
            
            
        end
        %New Normalized enrichment is added to the  existing ones
        fact_x = sqrt(trapz(x,Sx.^2)/(x(end)-x(1))^2);
        fact_t = sqrt(trapz(t,St.*St)/(t(end)-t(1))^2);
        fact_k = sqrt(trapz(k,Sk.*Sk)/(k(end)-k(1))^2);
        
        fact_xtk = (fact_x*fact_t*fact_k)^(1/3);
        X = [X fact_xtk*Sx/fact_x];
        T = [T fact_xtk*St/fact_t];
        K = [K fact_xtk*Sk/fact_k];
        
        
        %simplified stopping criterion (Eq. (2.26))
        E = sqrt(trapz(x,Sx.^2)*trapz(t,St.^2)*trapz(k,Sk.^2)) / sqrt(trapz(x,X(:,1).^2)*trapz(t,T(:,1).^2)*trapz(k,K(:,1).^2));
        if(E<epsilon_tilde), break; end;
        
    end
end
