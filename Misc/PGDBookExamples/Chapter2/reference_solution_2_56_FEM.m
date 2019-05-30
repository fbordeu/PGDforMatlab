function UU = reference_solution_2_56_FEM(x,y,GX,GY,q,Fx,Fy)
%Finite element implementation of the solution of problem 2.56
%
%Copyright (c) 2013, Francisco Chinesta (Ecole Centrale de Nantes), Roland Keunings (Universite catholique de Louvain), Adrien Leygue(CNRS)
%Author: Adrien Leygue.
%All rights reserved.
%See License file for more details


    %definition of mesh sizes
    hx = x(2)-x(1);
    hy = y(2)-y(1);
    
    %Mesh generation
    Domain = create_2Dmesh(x,y);
    
    %FEM matrix & RHS initialization
    K = sparse(Domain.Nnodes,Domain.Nnodes);
    B = zeros(Domain.Nnodes,1);
    
    %Elemental 2D Mass matrix
    Me = [  (hx*hy)/9  (hx*hy)/18 (hx*hy)/36 (hx*hy)/18; ...
        (hx*hy)/18 (hx*hy)/9  (hx*hy)/18 (hx*hy)/36; ...
        (hx*hy)/36 (hx*hy)/18 (hx*hy)/9  (hx*hy)/18; ...
        (hx*hy)/18 (hx*hy)/36 (hx*hy)/18 (hx*hy)/9];
    
    %Elemental 2D Stiffness matrix
    Ke = [   hx/(3*hy)+hy/(3*hx)    hx/(6*hy)-hy/(3*hx)  -hx/(6*hy)-hy/(6*hx)    hy/(6*hx)-hx/(3*hy); ...
        hx/(6*hy)-hy/(3*hx)    hx/(3*hy)+hy/(3*hx)   hy/(6*hx)-hx/(3*hy)   -hx/(6*hy)-hy/(6*hx); ...
        -hx/(6*hy)-hy/(6*hx)    hy/(6*hx)-hx/(3*hy)   hx/(3*hy)+hy/(3*hx)    hx/(6*hy)-hy/(3*hx); ...
        hy/(6*hx)-hx/(3*hy)   -hx/(6*hy)-hy/(6*hx)   hx/(6*hy)-hy/(3*hx)    hx/(3*hy)+hy/(3*hx)];
    
    %Elemental 1D mass matrix
    Me_bnd_x = [ hx/3 hx/6 ; hx/6, hx/3];
    
    %reconstruction of the 2D source term
    F = Fy*Fx';
    %reshape in the global node vector numbering
    F = F(:);
    
    %reconstruction of the 2D function providing the Dirichlet BC
    G = GY*GX';
    %reshape in the global node vector numbering
    G = G(:);
    
    %Assembly loop
    for i=1:Domain.Nelements
        %Gobal numbering of the nodes of the 2D element
        connectivity = Domain.elements(i,:);
        %Stiffness matrix
        K(connectivity,connectivity) = K(connectivity,connectivity) - Ke;
        %Source term
        B(connectivity) = B(connectivity) + Me*F(connectivity,1);
    end
    
    %Assembly of the Neumann boundary contribution
    for i = 1:size(Domain.segments_up,1)
        %Global numbering of the nodes of the 1D element
        connectivity = Domain.segments_up(i,:);
        %Neumann BC
        B(connectivity) = B(connectivity) - Me_bnd_x*(q*ones(2,1));
    end
    
    %Identification of the nodes with a Dirichlet BC
    Dirichlet = unique([Domain.nodes_left(:);Domain.nodes_right(:);Domain.nodes_down(:)]);
    %Simple imposition of the BC
    K(Dirichlet,:) = 0;
    K(Dirichlet,Dirichlet) = speye(numel(Dirichlet));
    B(Dirichlet) = G(Dirichlet);
    
    %System solution
    U = K\B;
    
    %reshape for visualization. Uses the fact that the nodes are generated
    %through the meshgrid command
    UU = reshape(U,[numel(y) numel(x)]);
end


function domain = create_2Dmesh(x,y)
    %size of the arguments
    Nx = numel(x);
    Ny = numel(y);
    
    %generate the nodes and their numbering
    [XX,YY] = meshgrid(x,y);
    domain.X = [XX(:) YY(:)];
    indices = reshape(1:(Nx*Ny),[Ny Nx]);
    
    %generate the elements (counterclockwise)
    nd1 = indices(1:end-1,1:end-1);
    nd2 = indices(1:end-1,2:end);
    nd3 = indices(2:end,2:end);
    nd4 = indices(2:end,1:end-1);
    domain.elements = [nd1(:) nd2(:) nd3(:) nd4(:)];
    
    %identify nodes on the different boundaries
    domain.nodes_down = indices(1,:)';
    domain.nodes_up = indices(end,:)';
    domain.nodes_left = indices(:,1);
    domain.nodes_right = indices(:,end);
    
    %identify segments of the different boundaries
    %!! segments have no prefered orientation !!
    domain.segments_down = [domain.nodes_down(1:end-1) domain.nodes_down(2:end)];
    domain.segments_up = [domain.nodes_up(1:end-1) domain.nodes_up(2:end)];
    domain.segments_left = [domain.nodes_left(1:end-1) domain.nodes_left(2:end)];
    domain.segments_right = [domain.nodes_right(1:end-1) domain.nodes_right(2:end)];
    
    domain.Nnodes = Nx*Ny;
    domain.Nelements = (Nx-1)*(Ny-1);
    
end
