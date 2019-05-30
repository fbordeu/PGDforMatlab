function result = sepdot(FF,GG)
%SEPDOT calculate the scalar product FF*GG : i.e. the sum of the product of
% the recosntruction of FF and GG
%
%   Note: there is no approximation in this operation
%
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
%

% is allNterms is bigger than dNterms a recursion is used to keep memory
% foot print low.

dNterms = 10000;

[Ndims,~,allNterms] = multiValidateFF(FF,GG);


if all(allNterms<=dNterms)
accum = FF{1}'*GG{1};
for d = 2:Ndims
    accum = accum.*(FF{d}'*GG{d});
end
result = sum(sort(accum(:)));
elseif allNterms(1)>dNterms
    limit = floor(allNterms(1)/2);
    tmp = cellfun(@(x) x(:,1:limit),FF,'uniformoutput',false);
    result = sepdot(tmp,GG);
    tmp = cellfun(@(x) x(:,(limit+1):allNterms(1)),FF,'uniformoutput',false);
    result = result+sepdot(tmp,GG);
elseif allNterms(2)>dNterms
    limit = floor(allNterms(2)/2);
    tmp = cellfun(@(x) x(:,1:limit),GG,'uniformoutput',false);
    result = sepdot(FF,tmp);
    tmp = cellfun(@(x) x(:,(limit+1):allNterms(2)),GG,'uniformoutput',false);
    result = result+sepdot(FF,tmp);
end

end