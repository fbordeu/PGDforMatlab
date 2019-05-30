function [Ndims,Ndofs,allNterms] = multiValidateFF(varargin)
%[Ndims,Ndofs,allNterms] = multiValidateFF(A,B,...)
%Checks if all inputs are valid separated representations over the same
%space (i.e. checks is all inputs have the same number of dimensions and the same 
%number dofs for each dimension) 
%
%Return values:
% Ndims: number of dimensions
% Ndofs: number of degrees of freedom per dimension (column vector)
% allNterms: number of terms in the separated representation of each input (line vector)
%
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Adrien Leygue (Adrien.Leygue@ec-nantes.fr)
%

if isempty(varargin)
    error('not enough input');
end
%allocate storage
NFF = numel(varargin);
allNdims = zeros(1,NFF);
allNterms = zeros(1,NFF);
allNdofs = cell(1,NFF);

%check individual inputs
for i=1:NFF
    [allNdims(i),allNdofs{i},allNterms(i)] = validateFF(varargin{i});
end

%check if the number of dimensions are the same between inputs
assert(all(allNdims==allNdims(1)),'Inconsistent number of dimensions between inputs');
Ndims = allNdims(1);

%checks if the number of dofs of each dimension is the same between inputs
allNdofs = horzcat(allNdofs{:});
tmp = bsxfun(@eq,allNdofs,allNdofs(:,1));
assert(all(tmp(:)),'Inconsistent number of dofs between inputs');
Ndofs = allNdofs(:,1);

end