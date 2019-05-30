function [s,n,d] = sepexppade( arg1,varargin )
%SEPEXP exponencial of a separated field using a Pade expansion
%%        Only valid for values between [-Pi Pi]
%%
%%  s = sepcos( FF, options,... )
%%
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
terms = 10;

for k = 1:length(varargin)
    if strcmpi(varargin{k}, 'order')
        terms = ceil(varargin{k+1}/2);
        varargin = varargin([1:(k-1) (k+2):end]);
        break
    end
end


p = terms;
q = terms;

for i = 0:p
    %%h(p+1-i) = (factorial(p)/factorial(p-i))/(factorial(i)*factorial(p+q)/factorial(p+q-i));
    h(p+1-i) = pfactorial(p-i+1,p)/(factorial(i)*pfactorial(p+q-i+1,p+q));
end

for j = 0:q
    %%k(q+1-j) = (-1)^(j)*(factorial(q)/factorial(q-j))/(factorial(j)*factorial(p+q)/factorial(p+q-j));
    k(q+1-j) = (-1)^(j)*pfactorial(p-j+1,p)/(factorial(j)*pfactorial(p+q-j+1,p+q));
end
%h = [1/840 2/105 1/7 4/7 1];
%k = [-1/210 1/14 -3/7 1];

n = seppolyval( h, arg1, varargin{:});
d = seppolyval( k, arg1, varargin{:});


s = quotient(n,d);


end

function res =pfactorial(start,endd)

res =1;
for i = start:endd
    res = res *i;
end

end 

