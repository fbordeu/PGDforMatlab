function err = Compute_quadratic_error(FF1,FF2,N_NT)
%This function computes the quadratic error between FF1 and FF2 using the
%integration matrices in N_NT
%Outputs: 
%       err :  multidimensional integral of (FF1-FF2)^2 computed using the mass matrices found in N_NT 
%Inputs:
%       FF1: separated representation of a multidimensional function
%       FF1: separated representation of a multidimensional function
%       N_NT: cell array of matrices containing the mass matrix along each
%       dimension
%
%Copyright (c) 2013, Francisco Chinesta (Ecole Centrale de Nantes), Roland Keunings (Universite catholique de Louvain), Adrien Leygue(CNRS)
%Author: Adrien Leygue.
%All rights reserved.
%See License file for more details


%extract relevant scalar values
%dimensionality
dim = numel(N_NT);
%number of terms in FF1
nf1 = size(FF1{1},2);
%number of terms in FF2
nf2 = size(FF2{1},2);

%Initializes a container for computing the integrals of all the cross products 
%of the different terms the particular initialization will yield the
%integral of the square of the difference between FF1 and FF2
tmp = [ones(1, nf1) -ones(1,nf2)];
prod_aux = tmp'*tmp;

%Compute the multidimensional integral of all the terms in the separated
%representation of (FF1-FF2)^2
for dd = 1:dim
    prod_aux(1:nf1,1:nf1) = prod_aux(1:nf1,1:nf1) .* (FF1{dd}' *  N_NT{dd} * FF1{dd});
    prod_aux((nf1+1):(nf1+nf2),(nf1+1):(nf1+nf2)) = prod_aux((nf1+1):(nf1+nf2),(nf1+1):(nf1+nf2)) .* (FF2{dd}' *  N_NT{dd} * FF2{dd});
    prod_aux(1:nf1,(nf1+1):(nf1+nf2)) = prod_aux(1:nf1,(nf1+1):(nf1+nf2)) .* (FF1{dd}' *  N_NT{dd} * FF2{dd});
end
prod_aux((nf1+1):(nf1+nf2),1:nf1) = (prod_aux(1:nf1,(nf1+1):(nf1+nf2)))';
%Sum all contributions
err = (abs(sum(sort(prod_aux(:)))));
end