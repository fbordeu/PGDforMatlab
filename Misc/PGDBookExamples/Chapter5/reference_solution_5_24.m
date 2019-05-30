function [U] = reference_solution_5_24(t,x,k,N)
    %This function computes the expression (5.24) of the book on the 2D grid
    %generated by the 1D grids t and x for the provided value of k.
    %A total of N term of the solution is computed
    %
    %Copyright (c) 2013, Francisco Chinesta (Ecole Centrale de Nantes), Roland Keunings (Universite catholique de Louvain), Adrien Leygue(CNRS)
    %Author: Adrien Leygue.
    %All rights reserved.
    %See License file for more details
    
    
    %auxilliary variable for convenience
    L = x(end);
    
    %indices in the summation
    n = 1:2:(2*(N-1));
    
    %Computation of the 1D functions in the solution
    X =bsxfun(@(x,n) sin(n*pi*x/L),x(:),n);
    T =bsxfun(@(t,n) exp(-k*n^2*pi^2*t/L^2),t(:),n);
    
    %scalar factor for each term
    A = -4*L^3/(k*pi^3)*n.^(-3);
    
    %global solution
    U = bsxfun(@times,X,A)*T' +  bsxfun(@times,(x(:).*(L-x(:)))/2/k,ones(1,numel(t)));
end