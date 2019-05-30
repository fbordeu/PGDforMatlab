function [alphas, varargout] = getalphas(FF, P )
% GETALPHAS extract alphas for a given separated field
%
%  [alphas,[FFnormalized]] = getalphas( FF, [P=2] )
%
%  P is the kind of norm used for the calculation (default 2).
%  See the documentation of norm
%
% See applyalphas, norm.
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
%

%%Validation
[Ndims,Ndofs,Nterms] = validateFF(FF);

%%Options
if~exist('P') 
    P = 2;
end

%% Initialization
alphas = ones(1,Nterms);

%% calcul of alphas
for dim=1:Ndims
    alphas_tmp = sum(abs(FF{dim}).^P).^(1/P);
    
    % normalization
    if (nargout > 1)  
        FF{dim} = bsxfun(@rdivide, FF{dim}, alphas_tmp);
    end

    alphas = alphas.*alphas_tmp;
end

%% Output Generation
if (nargout > 1)
    varargout{1} = FF;
end

