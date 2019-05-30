function res = PopRecompactOptions(varargin)
% POPRECOMPACTOPTIONS to pop a recompact structure of options from the pile
% 
% function res = PopRecompactOptions([true|false])
%
%  return the the poped structure
%
%   if false is use as argument the pile will not be modified
%
% See also PushRecompactOptions, CleanRecompactOptions, recompact.
%
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
%
global ReCompactOptions_DoNotRemove


if numel(ReCompactOptions_DoNotRemove) ~= 0
    res = ReCompactOptions_DoNotRemove(end);
    
    if nargin == 1 && varargin{1} == false
        return 
    end
    
    %disp('removing')
    ReCompactOptions_DoNotRemove = ReCompactOptions_DoNotRemove(1:end-1);
    if numel(ReCompactOptions_DoNotRemove) == 0
        %disp('undefining1')    
        clear global ReCompactOptions_DoNotRemove;    
    end
    
else
    %fprintf(2,'Cant Pop Recompact Options\n')
    clear global ReCompactOptions_DoNotRemove;
    res = recompact;
end
