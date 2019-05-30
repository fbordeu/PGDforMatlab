function  Y  = sepeval( arg1, Index )
%SEPEVAL to evaluate a separated field in one point
%   
% Y  = sepeval( FF, Index )
%
%    Note: there is no approximation in this operation
%
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
%
[Ndims,Ndofs,Nterms] = validateFF(arg1);

for i = 1: Ndims
    assert(~any(Ndofs(i) < Index(:,i)),'ERROR Index out of range');
end

pre = ones(size(Index,1),Nterms);

for i = 1:Ndims
   pre = pre .* arg1{i}(Index(:,i),:);
end

Y = sum(pre,2);

end

