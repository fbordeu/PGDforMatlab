function result = sepnorm(FF)
%SEPNORM calculate the norm of FF : i.e. the square root of the sum of the 
%        reconstruction of FF
%
%  function result = sepnorm(FF)
%
%   Note: there is no approximation in this operation
%
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
%

% if allNterms is bigger than dNterms a recursion is used to keep memory
% foot print low.

dNterms = 10000;

[Ndims,~,Nterms] = validateFF(FF);

if Nterms<dNterms
    accum = FF{1}'*FF{1};
    for d = 2:Ndims
        accum = accum.*(FF{d}'*FF{d});
    end
    result = sum(sort(accum(:)));
else
    limit = floor(Nterms(1)/2);
    tmp1 = cellfun(@(x) x(:,1:limit),FF,'uniformoutput',false);
    tmp2 = cellfun(@(x) x(:,(limit+1):Nterms),FF,'uniformoutput',false);
    result = sepnorm(tmp1)^2 +sepnorm(tmp2)^2 + 2*sepdot(tmp1,tmp2);
end
result = sqrt(max([result 0]));

end