%FEM_SOLVE routine
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Adrien Leygue (Adrien.Leygue@ec-nantes.fr)
%
function result = FEM_SOLVE(A,b,Dirichlet_pos,Dirichlet_val,varargin)
    N = size(A,1);
    mask = true(N,1);
    mask(Dirichlet_pos) = false;
    [~,p] = sort(Dirichlet_pos);
    Dirichlet_val = Dirichlet_val(p,:);
    result = zeros(N,size(b,2));
    result(~mask,:) = Dirichlet_val;
    if isempty(Dirichlet_pos)
        result(mask,:) = A(mask,mask) \ b(mask,:);
    else
        result(mask,:) = A(mask,mask) \ (b(mask,:) - A(mask,~mask)*Dirichlet_val);
    end
    
end
