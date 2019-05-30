function  s = sepcos( arg1, varargin )
% SEPCOS Cosinus of a separated field using a Taylor expansion
%        Only valid for values between [-2Pi 2Pi]
%  
%  s = sepcos( FF, options,... )
%
%  options are options for the internal recompact.
%  and/or order of the expation ('order',20)
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
%power 2 and rest
% we work for each factor
prod2 = sepprod(arg1,arg1,varargin{:});
prod = sepprod(0.5,prod2,varargin{:});

s = sepsum(1,sepprod(-1,prod,varargin{:}),varargin{:});

for k = 2:terms;
   p = 1/(2*k*(2*k-1));
   disp(['working in power ' num2str(k)])
   prod = sepprod(sepprod(p,prod2,varargin{:}),prod,varargin{:});
   s = sepsum(s,sepprod((-1)^k,prod,varargin{:}),varargin{:});
end

end

