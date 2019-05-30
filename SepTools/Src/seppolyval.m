function s = seppolyval( p, arg1, varargin)
%SEPPOLYVAL Evaluate polynomial of a separated field with internal recompaction
%    Y = seppolyval(P,X, options) returns the value of a polynomial P 
%    evaluated at X. P is a vector of length N+1 whose elements are the 
%    coefficients of the polynomial in descending powers.
% 
%        Y = P(1)*X^N + P(2)*X^(N-1) + ... + P(N)*X + P(N+1)
%
%  options are options for the internal recompact.
%
% See also PushRecompactOptions, PopRecompactOptions, sepsum, sepprod, recompact.
% 
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
%

[Ndims,Ndofs,Nterms] = validateFF(arg1);

s = cell(Ndims,1);


% initialization 
for d = 1:Ndims
  s{d} = zeros(Ndofs(d),0 );
end

% power 0, addition of the constant
s = sepsum(p(end),s,varargin{:});

% power 1
if numel(p) == 1
   return
end

if (p(end-1) ~= 0)
   s = sepsum(s,sepprod(p(end-1),arg1,varargin{:}),varargin{:});
end

%power 2 and the rest
% we work for each factor
prod = arg1;
for x = numel(p)-2:-1:1
  deg = numel(p)-x;
  disp(['working in power' num2str(deg)])
  prod = sepprod(prod,arg1,varargin{:});
  if (p(x) ~= 0)
     s = sepsum(s,sepprod(p(x),prod,varargin{:}),varargin{:});
  end
end

end



