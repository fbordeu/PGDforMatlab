function [ s ] = sepfracpow(arg1, alpha, varargin)
%sepfracpow fractional power of a separated
%   Using a power serie.
%   WARNING this serie has a radius of convergence of 1
%
%  [ s ] = sepfracpow(arg1,alpha, options )
%
%  s = agr1^(alpha) around x0
%
%  NOTE: The field MUST be positive 
%   FF always > 0  
%
%  options are options for the internal recompact.
%  and/or 
%  the order of the expation ('order',30)
%  the starting point of the expation ('around',1)
%
% 
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
%
validateFF(arg1);

dims1 =   numel(arg1);

order = 30;
x0=1;

for k = 1:length(varargin)
    if strcmpi(varargin{k}, 'order')
        order = varargin{k+1};
        varargin = varargin([1:(k-1) (k+2):end]);
        break
    end
end

for k = 1:length(varargin)
    if strcmpi(varargin{k}, 'around')
        x0 = varargin{k+1};
        varargin = varargin([1:(k-1) (k+2):end]);
        break
    end
end

s = cell(dims1,1);
% power 0
% initialization and addition of the constant
for d = 1:dims1
    s{d} = ones(size(arg1{d},1),1 );
end
s = sepprod(s,x0^alpha);


p = 1;
prod = 1;
for k = 1:order
    p = p*(alpha-k+1)/k;
    f = p*x0^(alpha-k);
    % in the case we give a integer as alpha
    if f == 0
        break
    end
    prod = sepprod(prod,sepsum(-1*x0,arg1,varargin{:}),varargin{:});
    s = sepsum(s,sepprod(f,prod,varargin{:}),varargin{:});
end


end
