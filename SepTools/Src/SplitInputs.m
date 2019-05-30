function [args1, args2]= SplitInputs(input)
% SPLITINPUST split the input cell vector into 2 cells vector
%
%    function [args1, args2]= SplitInputs(input)
%
%  args1 contains only the firsts cells and doubles
%
%  args2 contains the rest of the element of the input cells
%
%  See also sepsum.
%
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
%
   mask = cellfun(@(arg)((iscell(arg) || isscalar(arg)) && ~isstruct(arg)),input);
   l = find(~mask, 1);
   l = [l numel(mask)+1];
   args1 = input(mask(1:l(1)-1));
   args2 = input(l(1):end);
    
end