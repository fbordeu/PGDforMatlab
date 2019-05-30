function [s,n,d] = sepmexppade( arg1,varargin )
%SEPEXPMXPADE exponencial of -x^2 of a separated field using a Pade expansion
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


%v = ceil(log2(sepmax(sepprod(arg1,-1)) )/2);
%v = 2%
%arg1 = sepfrac(arg1,2^v, varargin{:});
pade.pn = [-4.257410245994015e-09  0.023862952603231  0.488945955642089  2.504590660104117];
pade.pd = [0.076148789064821 -0.105661545714008 0.855347930936226  -2.000676834671535 2.505024441014189];

%pade.pn = [0.032642005260451 0.120644981845740 0.263221818748067 0.326014934213662 0.183051101967515];
%pade.pd = [0.009606642854931 -0.110975942311374 0.487515479452915 -1 0.820377474304616];
%pade.min = -1.500000000000000;
%pade.max = 1.500000000000000;
    
s = seppadeval(pade, arg1, varargin{:});

%for i = 1:v;
%    s = sepprod(s,s);
%end 

end


