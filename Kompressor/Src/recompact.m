function res = recompact(datain, varargin)
%
%Funciton to recompact a separated solution
%
%recompact(filename, options,...)
%
%
%  datar = recompact('import_cas_30_0_parametres_500_modes.pxdmf', options,... )
%       using first argument as a filename
%
%  datar = recompact(data, options,... )
%       using first argument as a vectors of XdmfGrid
%
%  datar = recompact(data, options,... )
%       using first argument as a pxdmf data
%
%  fieldr = recompact(field, options,... )
%       using first argument as a cell array with dims matrices
%
%
%       'reweight_modes'    [true]|false
%       'fp_tol'            [1E-8]
%       'res_reduc'         [1e-8]
%       'max_added_modes'   [Number of modes in the solution]
%       'fp_max_iter'       [200]
%       'improve_modes'     [true]|false
%       'improve_modes_max' [10% of max_added_modes]
%       'usev6' to use the V6 version for compresion
%       'useTCPIP' to use the TCPIP connection to the compresion (experimental)
%       'save' to automaticaly save the solution the file to filenameR.pxdmf
%       'lastImproveModesLoop' [1] number of iteration of the improvement
%                                  loop on the last enrichement.
%       'sep_comp'         [true]|false to recompact each composant individualy
%       'verbose'          [true]|false
%       'donttouch'        {''} a cell of fields names to do not do the recompaction
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
%

%the time bomb
if all([ (clock >= [2019 12 0 0 0 0])  (builtin('clock') >= [2019 12 0 0 0 0])] )
  error('License Expired!!!');
end

opt.improve_modes     = 1;
opt.improve_modes_max = 0; % if zero is the the max(10, 10% of  the number of modes)
opt.useV6 = 0;
opt.reweight_modes    = 1;
opt.fp_tol            = 1E-8;
opt.res_reduc         = 1e-8;
opt.max_added_modes   = 0; % if zero the the max added modes is the number of modes in the original field
opt.fp_max_iter       = 200;
opt.verbose = true;
opt.lastImproveModesLoop = 1;
opt.sep_comp = true;
opt.recompact = true;

if (nargin == 0)
    if  (nargout == 0)
        help recompact
        return
    else
        res = opt;
        return ;
    end
end
opt = PopRecompactOptions(false);

opt.save = 0;
opt.recompile = false;
opt.debug = false;

opt.useTCPIP = 0;
opt.host = 'pc-lmm29';
opt.port = '5556';
opt.donttouch = {};



%% option parsing %%%%%%%%%%%%%%%%%%%%%%
for k = 1:2:length(varargin)
    if strcmpi(varargin{k}, 'sep_comp')
        opt.sep_comp = varargin{k+1};
        continue
    end
    if strcmpi(varargin{k}, 'reweight_modes')
        opt.reweight_modes = varargin{k+1};
        continue
    end
    if strcmpi(varargin{k}, 'fp_tol')
        opt.fp_tol = varargin{k+1};
        continue
    end
    if strcmpi(varargin{k}, 'res_reduc')
        opt.res_reduc = varargin{k+1};
        continue
    end
    if strcmpi(varargin{k}, 'max_added_modes')
        opt.max_added_modes = varargin{k+1};
        continue
    end
    if strcmpi(varargin{k}, 'fp_max_iter')
        opt.fp_max_iter = varargin{k+1};
        continue
    end
    if strcmpi(varargin{k}, 'improve_modes')
        opt.improve_modes = varargin{k+1};
        continue
    end
    if strcmpi(varargin{k}, 'improve_modes_max')
        opt.improve_modes_max = varargin{k+1};
        continue
    end
    if strcmpi(varargin{k}, 'useV6')
        opt.useV6 = varargin{k+1};
        continue
    end
    if strcmpi(varargin{k}, 'useTCPIP')
        opt.useTCPIP = varargin{k+1};
        continue
    end
    if strcmpi(varargin{k}, 'host')
        opt.host = varargin{k+1};
        continue
    end
    if strcmpi(varargin{k}, 'port')
        opt.port = varargin{k+1};
        continue
    end
    if strcmpi(varargin{k}, 'save')
        opt.save = varargin{k+1};
        continue
    end
    if strcmpi(varargin{k}, 'verbose')
        opt.verbose = varargin{k+1};
        continue
    end
    if strcmpi(varargin{k}, 'recompile')
        opt.recompile = varargin{k+1};
        continue
    end
    if strcmpi(varargin{k}, 'lastImproveModesLoop')
        opt.lastImproveModesLoop = varargin{k+1};
        continue
    end
    if strcmpi(varargin{k}, 'recompact')
        opt.recompact = varargin{k+1};
        continue
    end
    if strcmpi(varargin{k}, 'donttouch')
        opt.donttouch = varargin{k+1};
        continue
    end
    if strcmpi(varargin{k}, 'debug')
        opt.debug = varargin{k+1};
        if opt.debug 
            opt.verbose = true;
        end
        continue
    end
    error('Recompact:main:optnotfound', ['ERROR unknown option "' varargin{k} '"'])
end
%% to skip recompact
if opt.recompact == false
   res = datain;
   return
end
%% TCP IP connection %%%%%%%%%%%%%%%%%%%%%%%%%%
if opt.useTCPIP
    if(exist('pnet','file') == 0  || opt.recompile)
        if(opt.verbose); disp('need to compile the pnet.cpp file'); end
        type = computer;
        cpw = pwd();
        libpath = mfilename('fullpath');
        cd(libpath(1:find(libpath==filesep(),1,'last')))
        if strcmpi(type(1:4), 'GLNX')
            mex pnet.c
        else
            if strcmpi(type(1:3), 'MAC')
                %mex LDFLAGS='-bundle -framework accelerate $matlabroot/bin/maci64/libut.dylib'  'recompactmex.cpp'
                mex('CFLAGS=$CFLAGS -Wno-pointer-sign',['LDFLAGS=$LDFLAGS -bundle -framework accelerate ' matlabroot '/bin/maci64/libut.dylib'],'pnet.c')
            else
                disp('compiling in Windows')
            end
        end
        cd(cpw);
    end
    try
        opt.con=pnet('tcpconnect',opt.host,opt.port);
        if opt.con==-1; error 'Bad url or server down.....'; end
        if(opt.verbose); disp(['Connected to: ' opt.host]);end
    catch
        disp('ERROR in the connection')
        opt.useTCPIP = 0;
    end
end
%% C++ compilation %%%%%%%%%%%%%%%%%%%%%%%%%% 
if(opt.useV6==0 && opt.useTCPIP==0 )
    if(exist('recompactmex','file') == 0  || opt.recompile)
        if(opt.verbose); disp('need to compile the recompactmex.cpp file'); end
        type = computer;
        cpw = pwd();
        libpath = mfilename('fullpath');
        cd(libpath(1:find(libpath==filesep(),1,'last')))
        %to know witch blas is used call blaslib = fullfile(matlabroot, 'extern', 'lib', computer('arch'), 'microsoft', 'libmwblas.lib');
        if opt.debug
            mex -lmwblas -g -lut recompactmex.cpp cRecompactCore.cpp  PGD_Options.cpp                
        else
            if strcmpi(type(1:3), 'MAC')
              mex('-O',['CXXFLAGS=$CXXFLAGS -std=c++11'],['LDFLAGS=$LDFLAGS -bundle -lmwblas ' matlabroot '/bin/maci64/libut.dylib'],'recompactmex.cpp','cRecompactCore.cpp','PGD_Options.cpp')
            else
              mex -lmwblas -lut recompactmex.cpp cRecompactCore.cpp  PGD_Options.cpp
            end    
        end
  
%         if strcmpi(type(1:4), 'GLNX')
%             if opt.debug
%                 %mex -lgslcblas -g -lut recompactmex.cpp cRecompactCore.cpp  PGD_Options.cpp    
%             else
%                 %mex -lgslcblas -lut recompactmex.cpp cRecompactCore.cpp  PGD_Options.cpp
%             end
%         else
%             if strcmpi(type(1:3), 'MAC')
%                 %mex LDFLAGS='-bundle -framework accelerate /Applications/MATLAB_R2013a.app/bin/maci64/libut.dylib' recompactmex.cpp cRecompactCore.cpp
%                 mex('-O','-v',['CXXFLAGS=$CXXFLAGS -std=c++11'],['LDFLAGS=$LDFLAGS -bundle -framework accelerate ' matlabroot '/bin/maci64/libut.dylib'],'recompactmex.cpp','cRecompactCore.cpp','PGD_Options.cpp')
%             else
%                 disp('compiling in Windows')
%                 mex -lmwblas.lib -lut recompactmex.cpp cRecompactCore.cpp  PGD_Options.cpp
%             end
%         end
        cd(cpw);
    end
    if(exist('recompactmex','file') == 0)
        disp('compilation of mex failed using V6 for recompression')
        opt.useV6 = 1;
    end
end

%% type of the incomming data%%%%%%%%%%%%

switch class(datain)
    case { 'char' }
        res = recompact_file(datain, opt);
    case { 'struct' }
        if isfield(datain,'nodes')
            if(opt.verbose);  disp('using first argument as a pxdmf data');end
            res = recompact_data(datain,opt);
        else
            disp('Unknown struct type.')
        end
    case { 'cell' }
         if(numel(datain)==0)
            return;
         end
        validateFF(datain);
        if(opt.verbose); disp('using first argument as a field');end
        res = recompact_field(datain,opt);

    otherwise
        if isa( datain,'XdmfGridBase' )
            if datain.CheckIntegrity() == false ; 
                disp('Data Not well formed.')
                return;
            else
                res = RecompactXdmfGrid(datain,opt);
            end
        else
            disp('Unknown type for the data.')
            return
        end
end


if opt.useTCPIP
    pnet_putvar(opt.con,0);
    pnet(opt.con,'close');
end

if opt.save
    res.filename = [res.filename(1:end-6) 'R.pxdmf'];
    res.bin=1;
    res.precision = 'single';
    if(opt.verbose) ;disp(['saving file to : ' res.filename]);end;
    writepxdmf(res);
end

end

%% File compression
function res = recompact_file(filename,opt)
    if(opt.verbose); disp('using first argument as a filename'); end
    if(opt.verbose); disp(['Filename : ' filename]); end
    res = recompact_data(readpxdmf(filename),opt);
end

%% 
function res = recompact_data(data,opt)

if opt.sep_comp 
   datain = dataToScalar(data);
else
   datain = data;
end

dim = size(datain.nodes,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%   compaction of the node fields %%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:size(datain.nodes_fields,2)
    if(opt.verbose); disp([' compacting : ' datain.nodes_fields_names{i}]);end;
    clear solin
    solin = cell(dim,1);
    for j = 1:dim
        solin{j,1} = datain.nodes_fields{j,i}';
    end
    
    if ( ~any(strcmpi(datain.nodes_fields_names{i},opt.donttouch)) )
        solout = recompact_field(solin,opt);
        for j = 1:dim
            datain.nodes_fields{j,i} = solout{j}'  ;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   compaction of the cell fields %%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:size(datain.cell_fields,2)
    if(opt.verbose) ; disp([' compacting : ' datain.cell_fields_names{i}]);end;
    clear solin
    if ( ~any(strcmpi(datain.cell_fields_names{i},opt.donttouch)) )
        solin = cell(dim,1);
        for j = 1:dim
           solin{j} = datain.cell_fields{j,i}';
        end
        solout = recompact_field(solin,opt);
        for j = 1:dim
            datain.cell_fields{j,i} = solout{j}'  ;
        end
    end
end
res = datain;
end

function field = recompact_field(solin,opt)
if opt.useTCPIP == true
    field = quotient_field_TCPIP(solin,opt);
else
    field = recompact_field_Matlab(solin,opt);
end
if size(field{1},2) >= size(solin{1},2)
   if(opt.verbose); disp('You just kill a baby polar bear'); end
    disp('Warning Unable to recompact (returning the original solution)!!!!!')
    field = solin;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function field = quotient_field_TCPIP(FF,opt)

opt.max_added_modes

if(opt.max_added_modes ~=0)
    max_added_modes = opt.max_added_modes;
else
    max_added_modes = size(FF{1},2);
end
if (opt.improve_modes_max ~= 0)
    improve_modes_max = opt.improve_modes_max;
else
    improve_modes_max = max(10,round(size(FF{1},2)/10) );
end

options_t_be_sended = [opt.reweight_modes opt.fp_tol opt.res_reduc max_added_modes opt.fp_max_iter opt.improve_modes improve_modes_max opt.verbose opt.lastImproveModesLoop];

pnet_putvar(opt.con,2);

pnet_putvar(opt.con,options_t_be_sended);


NumberOfdims = numel(FF);
disp(['Sending number of dims' num2str(NumberOfdims )])
pnet_putvar(opt.con,NumberOfdims);
for i=1:numel(FF)
    %pause(1)
    b = double(FF{i});
    pnet_putvar(opt.con,b);
end
disp('Waiting for the results ')

field= cell(NumberOfdims,1);
for i=1:numel(FF)
    field{i} =  pnet_getvar(opt.con);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GridIn = RecompactXdmfGrid(GridIn,opt)


dim = numel(GridIn);

%% nodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nfields = GridIn(1).GetNumberOfNodeFields();
for i = 1:nfields
    name = GridIn(1).GetNodeFieldName(i);
    if ( ~any(strcmpi(name,opt.donttouch)) )
        clear solin
        solin = cell(dim,1);
        for j = 1:dim
            solin{j,1} =GridIn(j).GetNodeField(name);
        end
        solout = recompact_field(solin,opt);
        for j = 1:dim
             % we recover the number of the field 
             for l = 1: GridIn(j).GetNumberOfNodeFields() 
                 if strcmpi(GridIn(j).GetNodeFieldName(l),name); 
                     id = l; 
                 end; 
             end ;
             
             GridIn(j).nodeFields{id} = solout{j};
        end
    end
    
end

%% nodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nfields = GridIn(1).GetNumberOfElementFields();
for i = 1:nfields
    name = GridIn(1).GetElementFieldName(i);
    if ( ~any(strcmpi(name,opt.donttouch)) )
        clear solin
        solin = cell(dim,1);
        for j = 1:dim
            solin{j,1} =GridIn(j).GetElementField(name);
        end
        solout = recompact_field(solin,opt);
        for j = 1:dim
             % we recover the number of the field 
             for l = 1: GridIn(j).GetNumberOfElementFields() 
                 if strcmpi(GridIn(j).GetElementFieldName(l),name); 
                     id = l; 
                 end; 
             end ;
             
             GridIn(j).elementFields{id} = solout{j};
        end
    end
    
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function field = recompact_field_Matlab(solin,opt)

if(opt.max_added_modes ~=0)
    max_added_modes = opt.max_added_modes;
else
    max_added_modes = size(solin{1},2);
end
if (opt.improve_modes_max ~= 0)
    improve_modes_max = opt.improve_modes_max;
else
    improve_modes_max = max(10,round(size(solin{1},2)/10) );
end


if(opt.useV6==0 )
    for i =1:length(solin)
        solin{i} = double(solin{i});
    end
    field = recompactmex(solin,[opt.reweight_modes opt.fp_tol opt.res_reduc max_added_modes opt.fp_max_iter opt.improve_modes improve_modes_max opt.verbose opt.lastImproveModesLoop]);
else
    
    dims = length(solin);
    [problem, N_NT, options] = PGD_skeleton(dims,1);
    
    options.max_added_modes = max_added_modes;
    options.reweight_modes = opt.reweight_modes;
    options.improve_modes =  opt.improve_modes;
    options.fp_tol =  opt.fp_tol;
    options.improve_modes_max = improve_modes_max;
    options.fp_max_iter = opt.fp_max_iter;
    options.res_reduc =opt.res_reduc;
    options.verbose = opt.verbose;
    
    if opt.verbose
      disp('Using options : ')
      disp([' reweight_modes    ' num2str(options.reweight_modes ) ] );
      disp([' fp_tol            ' num2str(options.fp_tol ) ] );
      disp([' res_reduc         ' num2str(options.res_reduc ) ] );
      disp([' max_added_modes   ' num2str(options.max_added_modes ) ] );
      disp([' fp_max_iter       ' num2str(options.fp_max_iter ) ] );
      disp([' improve_modes     ' num2str(options.improve_modes ) ] );
      disp([' improve_modes_max ' num2str(options.improve_modes_max ) ] );
    end
    
    for i=1:dims
        problem.AA{i} = speye(size(solin{i},1));
        N_NT{i} =  speye(size(solin{i},1));
        problem.BB{i} = solin{i};
    end
    
    [sol2, err_reduc] = PGD_v6(problem,N_NT,options);
    
    for i=1:size(sol2.FF{1},2)
        sol2.FF{1}(:,i) = sol2.FF{1}(:,i)*sol2.alpha(i);
    end
    
    field = sol2.FF;

    display(['error reduction  : ' num2str(err_reduc(:)') ]);
    
end

end

%function printrecompactmex_cpp()
%fid = fopen('printrecompactmexcpp','w');
%fprintf(fid,'#include "mex.h" \n');
%fclose(fid);
%end

