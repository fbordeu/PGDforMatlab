function [ s ] = sepexp( arg1,varargin )
%SEPEXP exponencial of a separated field using a Taylor expansion
%        Only valid for values between [-Pi Pi]
%
%  s = sepcos( FF, options,... )
%
%  options are options for the internal recompact.
%  and/or order of the expation ('order',17)
%
% See also PushRecompactOptions, PopRecompactOptions, recompact.
%
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
%
terms = 17;

for k = 1:length(varargin)
    if strcmpi(varargin{k}, 'order')
        terms = varargin{k+1};
        varargin = varargin([1:(k-1) (k+2):end]);
        break
    end
end

% initialization 
% power 1

%power 2 and rest
% we work for each factor
s = sepsum(1,arg1,varargin{:});
prod = arg1;
for k = 2:terms;
  p = 1/k;
  deg = k;
  disp(['working in power' num2str(deg)])
  prod = sepprod(sepprod(p,arg1,varargin{:}),prod,varargin{:});
  s = sepsum(s,prod,varargin{:});
end


end

