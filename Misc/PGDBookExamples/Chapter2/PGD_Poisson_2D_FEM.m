function [X,Y] = PGD_Poisson_2D_FEM(x,y,Max_terms,Max_fp_iter,epsilon,epsilon_tilde,GX,GY,q,Fx,Fy)
    %PGD solution of problem 2.56, 2.57
    %For this example, we use Finite element technology
    %Outputs:
    %       X,Y :  computed PGD solution. The ith column of X (resp. Y) contains the nodal values of X_i (resp. Y_i)
    %              The first columns od X (resp. Y) contains GX (resp. GY).
    %Inputs:
    %       x,y: uniform 1D grids used along each dimension
    %       Max_terms: Maximum number of enrichments
    %       Max_fp_iter: Maximum number of fixed point iterations
    %       epsilon: termination criterion for the fixed point (Eq. (2.8))
    %       epsilon_tilde: termination criterion used for the enrichment process (Eq. (2.26))
    %       GX, GY: A priori known terms in the PGD solution to enforce
    %       non-homogeneous Dirichlet boundary conditions (Eq. (2.58))
    %       q: flux at the top boundary (Eq. (2.57), last line)
    %       Fx,Fy: Separated representation of the source term (Eq. (2.59))
    %
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
    
    %elemental stiffness matrices
    Kex = reshape([1.0./hx,-1.0./hx,-1.0./hx,1.0./hx],[2,2]);
    Key = reshape([1.0./hy,-1.0./hy,-1.0./hy,1.0./hy],[2,2]);
    
    %elemental mass matrices
    Mex = reshape([hx.*(1.0./3.0),hx.*(1.0./6.0),hx.*(1.0./6.0),hx.*(1.0./3.0)],[2,2]);
    Mey = reshape([hy.*(1.0./3.0),hy.*(1.0./6.0),hy.*(1.0./6.0),hy.*(1.0./3.0)],[2,2]);
    
    %Mass matrix along each dimension
    Mx = sparse(Nx,Nx);
    My = sparse(Ny,Ny);
    %Stiffness matrix along each dimension
    Kx = sparse(Nx,Nx);
    Ky = sparse(Ny,Ny);
    
    %Assembly of the FE matrices
    for el = 1:Nx-1
        Mx([el el+1],[el el+1]) = Mx([el el+1],[el el+1]) + Mex;
        Kx([el el+1],[el el+1]) = Kx([el el+1],[el el+1]) + Kex;
    end
    for el = 1:Ny-1
        My([el el+1],[el el+1]) = My([el el+1],[el el+1]) + Mey;
        Ky([el el+1],[el el+1]) = Ky([el el+1],[el el+1]) + Key;
    end
    
    %Account for the a priori known modes satisfying the non-homogeneous
    %Dirichlet BC.
    X = GX;
    Y = GY;
    %Number of a priori known terms
    N0 = size(X,2);
    
    %main enrichment loop
    for term=(N0+1):(Max_terms+N0)
        %initialization of the fixed point loop
        Sx = randn(Nx,1);
        Sy = randn(Ny,1);
        
        %Homogeneous Dirichlet Boundary conditions for Sx in 0 and Lx
        %Homogeneous Dirichlet Boundary conditions for Sy in 0 only
        Sx(1) = 0;
        Sx(end)=0;
        Sy(1) = 0;
        
        %fixed point iterations
        for iter=1:Max_fp_iter
            %Store the old values of SX & Sy for later comparison
            Sx_old = Sx;
            Sy_old = Sy;
            
            %Solve for Sx
            %construction of the boundary value problem along x
            %LHS coefficients
            alpha_x = Sy'*My*Sy;
            beta_x  = -Sy'*Ky*Sy;
            
            %Construction of the RHS
            %Source term coefficient
            ksi_x_i = Sy'*My*Fy;
            %Neumann BC coefficient
            mu_x = -q*Sy(end);
            
            RHS = Mx*( (Fx*ksi_x_i') + mu_x*ones(Nx,1));
            %In case this is not the first enrichment, previously computed
            %terms are added to the RHS
            if (term>1)
                %RHS coefficients
                gamma_x_i = Sy'*My*Y;
                delta_x_i = -Sy'*Ky*Y;
                
                RHS = RHS +(Kx*X)*gamma_x_i' - (Mx*X)*delta_x_i';
            end
            
            %construction of the FE boundary value problem
            A = -alpha_x*Kx + beta_x*Mx;
            %solution with homogeneous boundary conditions at both ends
            Sx(2:end-1) = A(2:end-1,2:end-1)\RHS(2:end-1);
            
            %Solve for Sy
            %construction of the boundary value problem along x
            %LHS coefficients
            alpha_y = Sx'*Mx*Sx;
            beta_y  = -Sx'*Kx*Sx;
            
            %Construction of the RHS
            %Source term coefficient
            ksi_y_i = Sx'*Mx*Fx;
            %Neumann BC coefficient
            mu_y = -q*(Sx'*Mx*ones(Nx,1));
            
            RHS = My*(Fy*ksi_y_i');
            RHS(end) = RHS(end) + mu_y;
            %In case this is not the first enrichment, previously computed
            %terms are added to the RHS
            if (term>1)
                %RHS coefficients
                gamma_y_i = Sx'*Mx*X;
                delta_y_i = -Sx'*Kx*X;
                
                RHS = RHS +(Ky*Y)*gamma_y_i' - (My*Y)*delta_y_i';
            end
            
            %construction of the FE boundary value problem
            A = -alpha_y*Ky + beta_y*My;
            %solution with homogeneous boundary conditions at first y node.
            Sy(2:end) = A(2:end,2:end)\RHS(2:end);
            
            %Norm of the difference between the 2 fixed point iterations  (Eq. (2.8))
            S_difference = sqrt((Sx'*Mx*Sx)*(Sy'*My*Sy) + (Sx_old'*Mx*Sx_old)*(Sy_old'*My*Sy_old) - 2*(Sx'*Mx*Sx_old)*(Sy'*My*Sy_old));
            %fixed point exit test
            if(S_difference < epsilon), break; end
            
            
        end
        %New normalized enrichment is added to the  existing ones
        fact_x = sqrt((Sx'*Mx*Sx)/(x(end)-x(1))^2);
        fact_y = sqrt((Sy'*My*Sy)/(y(end)-y(1))^2);
        fact_xy = sqrt(fact_x*fact_y);
        X = [X fact_xy*Sx/fact_x];
        Y = [Y fact_xy*Sy/fact_y];
        
        
        %simplified stopping criterion (Eq. (2.26))
        E = sqrt((Sx'*Mx*Sx)*(Sy'*My*Sy)) / sqrt((X(:,N0+1)'*Mx*X(:,N0+1))*(Y(:,N0+1)'*My*Y(:,N0+1)));
        
        if(E<epsilon_tilde), break; end;
        
    end
end
