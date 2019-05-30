function [X,Y] = PGD_Poisson_2D(x,y,Max_terms,Max_fp_iter,epsilon,epsilon_tilde,X0,Y0)
    %PGD solution of problem 2.1 with homogeneous Dirichlet boundary conditions
    %For this example, we use simple second order FD with trapezoidal
    %integration rule along each dimension
    %Outputs:
    %       X,Y :  computed PGD solution. The ith column of X (res. Y) contains the nodal values of X_i (resp. Y_i)
    %Inputs:
    %       x,y: uniform 1D grids used along each dimension
    %       Max_terms: Maximum number of enrichments
    %       Max_fp_iter: Maximum number of fixed point iterations
    %       epsilon: termination criterion for the fixed point (Eq. (2.8))
    %       epsilon_tilde: termination criterion used for the enrichment process (Eq. (2.26))
    %       X0,Y0: a priori known terms in the PGD solution. Optional inputs
    %
    %Copyright (c) 2013, Francisco Chinesta (Ecole Centrale de Nantes), Roland Keunings (Universite catholique de Louvain), Adrien Leygue(CNRS)
    %Author: Adrien Leygue.
    %All rights reserved.
    %See License file for more details
    
    
    
    %Mesh definition for each dimension
    Nx = numel(x);
    Ny = numel(y);
    %reshape to ensure that x & y are column vectors
    x = x(:);
    y = y(:);
    %Mesh size
    hx = (x(2)-x(1));
    hy = (y(2)-y(1));
    
    
    %Account for the presence of a priori known modes.
    if nargin==6
        %Matrix containing the nodal values of Xi(x)
        %The i^th column of X contains the nodal values of Xi(x);
        X = zeros(Nx,0);
        %Matrix containing the nodal values of Xi(x)
        %The i^th column of X contains the nodal values of Xi(x);
        Y = zeros(Ny,0);
    elseif nargin==8
        X = X0;
        Y = Y0;
    else
        error('You should provide both X & Y, not just X');
    end
    
    %Discrete operators along each dimension
    %Identity matrix long each dimension
    Ix = speye(Nx);
    Iy = speye(Ny);
    %Finite Differences second order differentiation matrices
    D2x = spdiags([ones(Nx,1)  -2*ones(Nx,1)  ones(Nx,1)],[-1 0 1],Nx,Nx)/hx^2;
    D2y = spdiags([ones(Ny,1)  -2*ones(Ny,1)  ones(Ny,1)],[-1 0 1],Ny,Ny)/hy^2;
    
    
    %Source term in separated form
    %f(x,y) = Fx(x)*FY(y)
    f = -1;
    Fx = f*ones(Nx,1);
    Fy = ones(Ny,1);
    
    %Number of a priori known terms
    N0 = size(X,2);
    
    %main enrichment loop
    for term=(N0+1):(Max_terms+N0)
        %initialization of the fixed point loop
        Sx = randn(Nx,1);
        Sy = randn(Ny,1);
        
        %homogeneous boundary conditions
        Sx(1) = 0;
        Sx(end)=0;
        Sy(1) = 0;
        Sy(end) = 0;
        
        %fixed point iterations
        for iter=1:Max_fp_iter
            %Store the old values of Sx & Sy for later comparison
            Sx_old = Sx;
            Sy_old = Sy;
            
            %Solve for Sx
            %construction of the boundary value problem along x
            %LHS coefficients
            alpha_x = trapz(y,Sy.^2);
            beta_x  = trapz(y,Sy.*(D2y*Sy));
            
            %Construction of the RHS
            ksi_x = trapz(y,Sy.*Fy);
            RHS = ksi_x*Fx;
            %In case this is not the first enrichment, previous terms are added
            %to the RHS
            if (term>1)
                %RHS coefficients
                gamma_x_i = trapz(y,bsxfun(@times,Sy,Y));
                delta_x_i = trapz(y,bsxfun(@times,Sy,D2y*Y));
                
                RHS = RHS -(D2x*X)*gamma_x_i' - X*delta_x_i';
            end
            
            %construction of the FD boundary value problem
            A = alpha_x*D2x + beta_x*Ix;
            %solution with homogeneous boundary conditions
            Sx(2:end-1) = A(2:end-1,2:end-1)\RHS(2:end-1);
            
            %Solve for Sy
            %construction of the boundary value problem along x
            %LHS coefficients
            alpha_y = trapz(x,Sx.^2);
            beta_y  = trapz(x,Sx.*(D2x*Sx));
            
            %Construction of the RHS
            ksi_y = trapz(x,Sx.*Fx);
            RHS = ksi_y*Fy;
            %In case this is not the first enrichment, previous terms are added
            %to the RHS
            if (term>1)
                %RHS coefficients
                gamma_y_i = trapz(x,bsxfun(@times,Sx,X));
                delta_y_i = trapz(x,bsxfun(@times,Sx,D2x*X));
                
                RHS = RHS -(D2y*Y)*gamma_y_i' - Y*delta_y_i';
            end
            
            %construction of the FD boundary value problem
            A = alpha_y*D2y + beta_y*Iy;
            %solution with homogeneous boundary conditions
            Sy(2:end-1) = A(2:end-1,2:end-1)\RHS(2:end-1);
            
            
            %Norm of the difference between the 2 fixed point iterations (Eq. (2.8))
            S_difference = sqrt(trapz(x,Sx.^2)*trapz(y,Sy.^2) + trapz(x,Sx_old.^2)*trapz(y,Sy_old.^2) - 2*trapz(x,Sx.*Sx_old)*trapz(y,Sy.*Sy_old));
            %fixed point exit test
            if(S_difference < epsilon), break; end
            
        end
        %Normalized new enrichment is added to the  existing ones.
        fact_x = sqrt(trapz(x,Sx.^2)/(x(end)-x(1))^2);
        fact_y = sqrt(trapz(y,Sy.^2)/(y(end)-y(1))^2);
        fact_xy = sqrt(fact_x*fact_y);
        X = [X fact_xy*Sx/fact_x];
        Y = [Y fact_xy*Sy/fact_y];
        
        %simplified stopping criterion (Eq. (2.26))
        E = sqrt(trapz(x,Sx.^2)*trapz(y,Sy.^2)) / sqrt(trapz(x,X(:,N0+1).^2)*trapz(y,Y(:,N0+1).^2));
        if(E<epsilon_tilde), break; end;
        
    end
end