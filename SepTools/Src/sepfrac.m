function sol = sepfrac( arg1, arg2, varargin )
%SEPFRAC Quotient of two separated fields with internal recompact
%   
% s = sepfrac( FF1, FF2, options )
%
% FF1 can be a scalar or a separated field
%
%  options are options for the internal recompact.
%
% See also PushRecompactOptions, PopRecompactOptions, recompact.
%
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
%

% scalar numerator
if ~iscell(arg1)
    
  % scalar denominator
  if ~iscell(arg2)
    sol =   arg1/arg2;
    return
  else
    val = arg1;
    arg1 = cell(size(arg2));
    for i= 1:numel(arg1)
        arg1{i} = ones(size(arg2{i},1),1);
    end
    arg1{1} = arg1{1}*val;
  end
else
   if ~iscell(arg2)
    sol= sepprod(arg1,1/arg2);
    return 
   end
end

multiValidateFF(arg1,arg2);

sol = quotient(arg1,arg2,varargin{:});
