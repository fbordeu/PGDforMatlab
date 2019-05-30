function PushRecompactOptions(data)
% PUSHRECOMPACTOPTIONS to push a recompact structure of options to the pile
%
% function  s = PushRecompactOptions( options )
%
%  if you use this function, be aware of global variables behavior in
%  matlab. Use the CleanRecompactOptions() to clean the pile.
%
% See also PopRecompactOptions, CleanRecompactOptions, recompact .
%
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
%

global ReCompactOptions_DoNotRemove
    
if numel(ReCompactOptions_DoNotRemove) ~= 0
    %disp('adding')
    ReCompactOptions_DoNotRemove(end+1) = data;
else
    %disp('defining')
    ReCompactOptions_DoNotRemove = data;
end

