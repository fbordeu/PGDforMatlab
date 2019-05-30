function setMatlabTools(varargin)
% Setup of MatlabTools 
%
% setMatlabTools(options,...)
%
% Options are:
%
% 'nosave' in this way the modified path will not be save for future sessions
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
%

opt.save = true;
for k = 1:nargin
    if strcmpi(varargin{k},'nosave')
        opt.save = false;
    end
end

% to to the SetMatlabtolls directory do the recurtion and then return to
% the courrent path
cpw = pwd();
libpath = mfilename('fullpath');
cd(libpath(1:find(libpath==filesep(),1,'last')))
    
disp('Setting up the path to all the Src directories (including subdirectories)')


content = dir();

for i = 1:numel(content)
    d = content(i);
    if d.name(1) == '.'
        continue
    end
    
    if d.isdir && exist([d.name filesep 'Src' filesep],'dir')
        p = [pwd() filesep d.name filesep 'Src' filesep];        
        fprintf(1,' ****  Adding %s to the path   ****\n', [d.name filesep 'Src' filesep]);
        addpath(p)
        disp(p)
        AddSubDirectoriesPath([pwd() filesep d.name filesep 'Src' filesep ]);
        
    end
end
%Manualy  added paths 

addManualPath([ pwd() filesep() 'Misc' filesep 'PXDMF' filesep]);

if (opt.save)
    savepath
    disp('Path saved for future sessions')
end
cd(cpw);

disp('**** Setup Done ****')
end

function AddSubDirectoriesPath(ptt)
% Helper function

content = dir(ptt);
for i = 1:numel(content)
    d = content(i);
    % to skeep files and directories containing '_' (includes_mpi)
    if d.name(1) == '.' || sum(d.name == '_') > 0
        continue
    end
    if d.isdir 
        p = [ ptt d.name];
        addpath(p)
        disp(p)
        AddSubDirectoriesPath([d.name filesep]);
    end
end

end

function addManualPath(p)
    addpath(p)
    disp(p)
end