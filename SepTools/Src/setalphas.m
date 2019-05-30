function FF = setalphas(FF,alphas, varargin )
% SETALPHAS inject alphas into a given separated field
%
%  [FF] = getalphas( FF, alphas, options, ...  )
%
%  by default will inject the alphas in the smallest dimension.
%  use option ('first',true) to inject in the first
%
% See setalphas, norm.
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
%

%%Validation
[Ndims,Ndofs,Nterms] = validateFF(FF);
assert(size(alphas,1) == 1,'alphas must be a row vector');
assert(size(alphas,2) == Nterms,'Inconsistent number of alphas and number of terms in the field');

%%Options
if length(varargin)

    for k = 1:length(varargin)
        if strcmpi(varargin{k}, 'first') && varargin{k+1}
            injectdim = 1;
            continue
        end
        [~,injectdim] = min(Ndofs);
    end
else
    [~,injectdim] = min(Ndofs);
end

%% alpha injection
FF{injectdim} = bsxfun(@times,FF{injectdim}, alphas);

