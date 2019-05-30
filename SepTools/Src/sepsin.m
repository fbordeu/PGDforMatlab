function  s = sepsin( arg1, varargin )
% SEPSIN Sinus of a separated field using a Taylor expansion.
%        Only valid for values between [-2Pi 2Pi]
%
%  s = sepsin( FF, options,... )
%
%  options are options for the internal recompact.
%  and/or order of the expation ('order',15)
%  
% See also PushRecompactOptions, PopRecompactOptions, recompact.
%
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
%
[Ndims,Ndofs,Nterms] = validateFF(arg1);

terms = 15;

for k = 1:length(varargin)
    if strcmpi(varargin{k}, 'order')
        terms = varargin{k+1};
        varargin = varargin([1:(k-1) (k+2):end]);
        break
    end
end

% initialization 
% power 1
s = arg1;

%power 3 and rest
% we work for each factor
prod = arg1;
prod2 = sepprod(arg1,arg1,varargin{:});
for k = 3:2:2*terms;
  p = (-1)^k/(k*(k-1));
  deg = k;
  disp([' Working in power ' num2str(deg)])
  prod = sepprod(sepprod(p,prod,varargin{:}),prod2,varargin{:}); 
  s = sepsum(s,prod,varargin{:});
end


end

