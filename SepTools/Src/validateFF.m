function [Ndims,Ndofs,Nterms] = validateFF(FF)
%[Ndims,Nukns,Nterms] = validateFF(FF)
%Checks if FF is a valid separated representation
%
%Return values:
% Ndims: number of dimensions
% Ndofs: number of degrees of freedom per dimension (column vector)
% Nterms: number of terms in the separated representation
%
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Adrien Leygue (Adrien.Leygue@ec-nantes.fr)
%

%do we have a non-empty cell column vector?
validateattributes(FF,{'cell'},{'column','nonempty'},'','FF');

Ndims = numel(FF);
Nterms = zeros(Ndims,1);
Ndofs = zeros(Ndims,1);

%is the content of each cell a nonempty 2d array of floating point numbers?
for i=1:Ndims
    validateattributes(FF{i},{'double','single'},{'2d'},'',['FF{' num2str(i) '}']);
    [Ndofs(i),Nterms(i)] = size(FF{i});
end

%do all the dimensions have the same number of terms
assert(all(Nterms==Nterms(1)),'Inconsistent number of terms between the dimensions');
Nterms = Nterms(1);
end