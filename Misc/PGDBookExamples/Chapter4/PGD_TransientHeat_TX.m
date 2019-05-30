function [X,T] = PGD_TransientHeat_TX(t,x,Max_terms,Max_fp_iter,epsilon,epsilon_tilde,Ft,Fx)
    %PGD solution of problem 4.1 with homogeneous Dirichlet boundary conditions
    %and homogeneous initial condition
    %For this example, we use simple second order FD with trapezoidal
    %integration rule along the space dimension and an implicit euler scheme
    %along the time dimension
    %Outputs:
    %       X,T :  computed PGD solution. The ith column of X (resp. T) contains the nodal values of X_i (resp. T_i)
    %Inputs:
    %       t,x: uniform 1D grids used along each dimension
    %       Max_terms: Maximum number of enrichments
    %       Max_fp_iter: Maximum number of fixed point iterations
    %       epsilon: termination criterion for the fixed point (Eq. (2.8))
    %       epsilon_tilde: termination criterion used for the enrichment process (Eq. (2.26))
    %       Ft,Fx: separated representation of the source term
    %
    %Copyright (c) 2013, Francisco Chinesta (Ecole Centrale de Nantes), Roland Keunings (Universite catholique de Louvain), Adrien Leygue(CNRS)
    %Author: Adrien Leygue.
    %All rights reserved.
    %See License file for more details
    
    
    
    %thermal diffusivity
    k=1;
    
    %Mesh definition for each dimension
    Nt = numel(t);
    Nx = numel(x);
    %reshape to ensure that x & t are column vectors
    t = t(:);
    x = x(:);
    %Mesh size
    dt = (t(2)-t(1));
    hx = (x(2)-x(1));
    
    %Discrete operators along each dimension
    %Finite difference matrices in time
    %Identity matrix
    It = speye(Nt);
    %time differentiation matrix
    D1t = spdiags([-ones(size(t)) ones(size(t))]/dt,[-1 0],Nt,Nt);
    D1t(1,:) = 0;
    
    %Finite difference matrices in x
    %Identity matrix
    Ix = speye(Nx);
    %Finite difference second order differentiation matrix
    D2x = spdiags([ones(Nx,1)  -2*ones(Nx,1)  ones(Nx,1)],[-1 0 1],Nx,Nx)/hx^2;
    
    X = zeros(Nx,0);
    T = zeros(Nt,0);
    
    
    %main enrichment loop
    for term=1:Max_terms
        %initialization of the fixed point loop
        St = randn(Nt,1);
        Sx = randn(Nx,1);
        
        %Satisfaction of the homogeneous Dirichlet Boundary conditions for the
        %enrichments
        St(1) = 0;
        Sx(1) = 0;
        Sx(end)=0;
        
        %fixed point iterations
        for iter=1:Max_fp_iter
            %Store the old values of Sx & St for later comparison
            St_old = St;
            Sx_old = Sx;
            
            %Solve for Sx
            %construction of the boundary value problem along x
            %LHS coefficients
            alpha_x = trapz(t,St.^2);
            beta_x  = trapz(t,St.*(D1t*St));
            
            %Source term coefficient
            ksi_x = Fx*trapz(t,bsxfun(@times,St,Ft))';
            
            %Construction of the RHS
            RHS = ksi_x;
            %In case this is not the first enrichment, previous terms are added
            %to the RHS
            if (term>1)
                %RHS coefficients
                gamma_x_i = trapz(t,bsxfun(@times,St,T));
                delta_x_i = trapz(t,bsxfun(@times,St,D1t*T));
                
                RHS = RHS + k*(D2x*X)*gamma_x_i' - X*delta_x_i';
            end
            
            %construction of the FD boundary value problem
            A = -alpha_x*k*D2x + beta_x*Ix;
            %solution with homogeneous boundary conditions
            Sx(2:end-1) = A(2:end-1,2:end-1)\RHS(2:end-1);
            
            %Solve for St
            %construction of the initial value problem along t
            %LHS coefficients
            alpha_t = trapz(x,Sx.^2);
            beta_t  = trapz(x,Sx.*(D2x*Sx));
            
            %Source term coefficient
            ksi_t = Ft*trapz(x,bsxfun(@times,Sx,Fx))';
            
            %Construction of the RHS
            RHS = ksi_t;
            %In case this is not the first enrichment, previous terms are added
            %to the RHS
            if (term>1)
                %RHS coefficients
                gamma_t_i = trapz(x,bsxfun(@times,Sx,X));
                delta_t_i = trapz(x,bsxfun(@times,Sx,D2x*X));
                
                RHS = RHS -(D1t*T)*gamma_t_i' + k*T*delta_t_i';
            end
            
            %construction of the initial value problem
            A = alpha_t*D1t - (beta_t*k)*It;
            %solution with homogeneous initial condition.
            St(2:end) = A(2:end,2:end)\RHS(2:end);
            
            
            %Norm of the difference between the 2 fixed point iterations  (Eq. (2.8))
            S_difference = sqrt(trapz(x,Sx.^2)*trapz(t,St.*St) + trapz(x,Sx_old.^2)*trapz(t,St_old.^2) - 2*trapz(x,Sx.*Sx_old)*trapz(t,St.*St_old));
            %fixed point exit test
            if(S_difference < epsilon), break; end
            
            
        end
        %New normalized enrichments are added to the  existing ones
        fact_x = sqrt(trapz(x,Sx.^2)/(x(end)-x(1))^2);
        fact_t = sqrt(trapz(St.*St)/(t(end)-t(1))^2);
        fact_xt = sqrt(fact_x*fact_t);
        X = [X fact_xt*Sx/fact_x];
        T = [T fact_xt*St/fact_t];
        
        
        %simplified stopping criterion (Eq. (2.26))
        E = sqrt(trapz(x,Sx.^2)*trapz(t,St.^2)) / sqrt(trapz(x,X(:,1).^2)*trapz(t,T(:,1).^2));
        
        if(E<epsilon_tilde), break; end;
        
    end
end
