function [Y,I] = sepmin( arg1, varargin )
%SEPMIN  min of a separated field
%
% [Y,I] = sepmin( FF, options )
%
%  NOTE: The field MUST be positive 
%   FF always > 0  
%
%  options are options for the internal recompact.
%  and/or 
%  Number of min searched in the field ('Nmax',1)
%  Maximal number of iteration for the power algo ('Mmax_iters',100)
%
% Y the minimum
% I the index
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
%
validateFF(arg1);


Y = sepmax( arg1,'Nmax',1,varargin{:});

% 
arg = sepsum(sepprod(-1, arg1),Y(1));

[~,I] = sepmax( arg, varargin{:} );

Y = sepeval(arg1,I);

end