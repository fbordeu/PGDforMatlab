function res = quotient(datain1,datain2, varargin)
%
%Funciton to quotient a separated solution
%
%
%  fieldr = quotient(nominatorfield, denominatorfield, options,... )
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
%       'sep_comp'         [true]|false to quotient each composant individualy
%       'verbose'          [true]|false
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

if (nargin == 0)
    if  (nargout == 0)
        help quotient
        return
    end
end
% this is experimental

opt.useV6             = false;
opt.useEasyQuotient   = false;
opt.reweight_modes    = 1;
opt.fp_tol            = 1E-8;
opt.res_reduc         = 1e-8;
opt.max_added_modes   = 0; % if zero the the max added modes is the number of modes in the original field
opt.fp_max_iter       = 200;
opt.improve_modes     = 1;
opt.improve_modes_max = 0; % if zero is the the max(10, 10% of  the number of modes)
opt.verbose           = true;
opt.recompile         = false;
opt.lastImproveModesLoop = 1;

opt.useTCPIP = false;
opt.host = 'pc-lmm29';
opt.port = '5556';

%% option parsing %%%%%%%%%%%%%%%%%%%%%%
for k = 1:2:length(varargin)
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
    if strcmpi(varargin{k}, 'useEasyQuotient')
        opt.useEasyQuotient = varargin{k+1};
        continue
    end
    
    if strcmpi(varargin{k}, 'useV6')
        opt.useV6 = varargin{k+1};
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
    error('quotient:main:optnotfound', ['ERROR unknown option "' varargin{k} '"'])
end

%% TCP IP connection %%%%%%%%%%%%%%%%%%%%%%%%%%
if opt.useTCPIP
    if(exist('pnet','file') == 0  || opt.recompile)
        if(opt.verbose); disp('need to compile the pnet.cpp file'); end
        type = computer; 
        cpw = pwd();
        libpath = mfilename('fullpath');
        cd(libpath(1:find(libpath=='/',1,'last')))
        if strcmpi(type(1:4), 'GLNX')
            mex pnet.c
        else
            if strcmpi(type(1:3), 'MAC')
                %mex LDFLAGS='-bundle -framework accelerate $matlabroot/bin/maci64/libut.dylib'  'quotientmex.cpp'
                mex('CFLAGS=$CFLAGS -Wno-pointer-sign',['LDFLAGS=$LDFLAGS -bundle -framework accelerate ' matlabroot '/bin/maci64/libut.dylib'],'pnet.c')
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
    if(exist('quotientmex','file') == 0  || opt.recompile)
        if(opt.verbose); disp('need to compile the quotientmex.cpp file'); end
        type = computer;
        cpw = pwd();
        libpath = mfilename('fullpath');
        cd(libpath(1:find(libpath==filesep(),1,'last')))
        if strcmpi(type(1:3), 'MAC')
            mex('-O',['CXXFLAGS=$CXXFLAGS -std=c++11'],['LDFLAGS=$LDFLAGS -bundle -lmwblas ' matlabroot '/bin/maci64/libut.dylib'],'quotientmex.cpp','cQuotientCore.cpp','PGD_Options.cpp')
        else
            mex -lmwblas -lut quotientmex.cpp cQuotientCore.cpp PGD_Options.cpp
        end     
%         if strcmpi(type(1:4), 'GLNX')
%             %mex -lgslcblas -lut -g quotientmex.cpp cQuotientCore.cpp PGD_Options.cpp
%             mex -lgslcblas -lut quotientmex.cpp cQuotientCore.cpp PGD_Options.cpp
%         else
%             if strcmpi(type(1:3), 'MAC')
%                 %mex LDFLAGS='-bundle -framework accelerate /Applications/MATLAB_R2013a.app/bin/maci64/libut.dylib' quotientmex.cpp cquotientCore.cpp
%                 mex('-O',['LDFLAGS=$LDFLAGS -bundle -framework accelerate ' matlabroot '/bin/maci64/libut.dylib'],'quotientmex.cpp','cQuotientCore.cpp','PGD_Options.cpp')
%             else
%                 disp('compiling in Windows')
%                 mex -lmwblas.lib -lut quotientmex.cpp cQuotientCore.cpp PGD_Options.cpp
%             end
%         end
        cd(cpw);
    end
    if(exist('quotientmex','file') == 0)
        disp('compilation of mex failed using V6 for recompression')
        opt.useV6 = 1;
    end
end

%% type of the incomming data%%%%%%%%%%%%

switch class(datain1)
    case { 'cell' }
        if(numel(datain1)==0 || numel(datain2)==0)
            return;
        end
        if(opt.verbose); disp('using first argument as a field');end
        multiValidateFF(datain1,datain2);
        res = quotient_field(datain1,datain2,opt);
    case { 'double'}
        factor = datain1;
        validateFF(datain2);
        datain1 = cellfun(@(arg1) (ones(size(arg1,1),1) ), datain2 ,'UniformOutput',false );
        datain1{1}  =  datain1{1}*factor;
        if(opt.verbose); disp('using first scalar argument as a constant field with one term');end
        res = quotient_field(datain1,datain2,opt);
    otherwise
        
    %    if isa( datain,'XdmfGridBase' )
    %        res = QuotientXdmfGrid(datain,opt);
    %    else
            error('Unknown type for the data.')
            return
    %    end
end

if opt.useTCPIP
    pnet_putvar(opt.con,0);
    pnet(opt.con,'close');
end

end


function field = quotient_field(solin1,solin2,opt)
if opt.useTCPIP == true
    field = quotient_field_TCPIP(solin1,solin2,opt);
else
    field = quotient_field_Matlab(solin1,solin2,opt);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function field = quotient_field_TCPIP(solin1,solin2,opt)

opt.max_added_modes

if(opt.max_added_modes ~=0)
    max_added_modes = opt.max_added_modes;
else
    max_added_modes = max(size(solin1{1},2),size(solin2{1},2));
end
if (opt.improve_modes_max ~= 0)
    improve_modes_max = opt.improve_modes_max;
else
    improve_modes_max = max(max(10,round(size(solin1{1},2)/10)),round(size(solin2{1},2)/10)) ;
end

options_t_be_sended = [opt.reweight_modes opt.fp_tol opt.res_reduc max_added_modes opt.fp_max_iter opt.improve_modes improve_modes_max opt.verbose opt.lastImproveModesLoop];

pnet_putvar(opt.con,2);

pnet_putvar(opt.con,options_t_be_sended);


NumberOfdims = numel(solin1);
disp(['Sending number of dims' num2str(NumberOfdims )])
pnet_putvar(opt.con,NumberOfdims);
for i=1:numel(solin1)
    %pause(1)
    b = double(solin1{i});
    pnet_putvar(opt.con,b);
end

NumberOfdims = numel(solin2);
disp(['Sending number of dims' num2str(NumberOfdims )])
pnet_putvar(opt.con,NumberOfdims);
for i=1:numel(solin2)
    %pause(1)
    b = double(solin1{i});
    pnet_putvar(opt.con,b);
end

disp('Waiting for the results ')

field= cell(NumberOfdims,1);
for i=1:numel(FF)
    field{i} =  pnet_getvar(opt.con);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function GridIn = QuotientXdmfGrid(GridIn,opt)
% 
% dim = numel(GridIn);
% 
% %% nodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nfields = GridIn(1).GetNumberOfNodeFields();
% for i = 1:nfields
%     name = GridIn(1).GetNodeFieldName(i);
%     if ( ~any(strcmpi(name,opt.donttouch)) )
%         clear solin
%         solin = cell(dim,1);
%         for j = 1:dim
%             solin{j,1} =GridIn(j).GetNodeField(name);
%         end
%         solout = quotient_field(solin,opt);
%         for j = 1:dim
%              % we recover the number of the field 
%              for l = 1: GridIn(j).GetNumberOfNodeFields() 
%                  if strcmpi(GridIn(j).GetNodeFieldName(l),name); 
%                      id = l; 
%                  end; 
%              end ;
%              
%              GridIn(j).nodeFields{id} = solout{j};
%         end
%     end
%     
% end
% 
% %% nodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nfields = GridIn(1).GetNumberOfElementFields();
% for i = 1:nfields
%     name = GridIn(1).GetElementFieldName(i);
%     if ( ~any(strcmpi(name,opt.donttouch)) )
%         clear solin
%         solin = cell(dim,1);
%         for j = 1:dim
%             solin{j,1} =GridIn(j).GetElementField(name);
%         end
%         solout = quotient_field(solin,opt);
%         for j = 1:dim
%              % we recover the number of the field 
%              for l = 1: GridIn(j).GetNumberOfElementFields() 
%                  if strcmpi(GridIn(j).GetElementFieldName(l),name); 
%                      id = l; 
%                  end; 
%              end ;
%              
%              GridIn(j).elementFields{id} = solout{j};
%         end
%     end
%     
% end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function field = quotient_field_Matlab(solin1,solin2,opt)

if(opt.max_added_modes ~=0)
    max_added_modes = opt.max_added_modes;
else
    max_added_modes = size(solin1{1},2);
end
if (opt.improve_modes_max ~= 0)
    improve_modes_max = opt.improve_modes_max;
else
    improve_modes_max = max(10,round(size(solin1{1},2)/10) );
end


if(opt.useV6==0 && opt.useEasyQuotient ==0 )
    
    
    for i =1:length(solin1)
        solin1{i} = double(solin1{i});
        solin2{i} = double(solin2{i});
    end
    field = quotientmex(solin1,solin2,[opt.reweight_modes opt.fp_tol opt.res_reduc max_added_modes opt.fp_max_iter opt.improve_modes improve_modes_max opt.verbose opt.lastImproveModesLoop]);
else
    
    
    if opt.useEasyQuotient 
       if opt.verbose
         disp('Using options : ')
         %disp([' reweight_modes    ' num2str(options.reweight_modes ) ] );
         disp([' fp_tol            ' num2str(opt.fp_tol ) ] );
         disp([' res_reduc         ' num2str(opt.res_reduc ) ] );
         disp([' max_added_modes   ' num2str(opt.max_added_modes ) ] );
         disp([' fp_max_iter       ' num2str(opt.fp_max_iter ) ] );
         disp([' improve_modes     ' num2str(opt.improve_modes ) ] );
         disp([' improve_modes_max ' num2str(opt.improve_modes_max ) ] );
       end
       
       field = quotientEasy(solin1,solin2,opt);

    else
       dims = length(solin1);
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
         %disp([' reweight_modes    ' num2str(options.reweight_modes ) ] );
         disp([' fp_tol            ' num2str(options.fp_tol ) ] );
         disp([' res_reduc         ' num2str(options.res_reduc ) ] );
         disp([' max_added_modes   ' num2str(options.max_added_modes ) ] );
         disp([' fp_max_iter       ' num2str(options.fp_max_iter ) ] );
         %disp([' improve_modes     ' num2str(options.improve_modes ) ] );
         %disp([' improve_modes_max ' num2str(options.improve_modes_max ) ] );
       end
       
       for i=1:dims
         for j=1:size(solin2{1},2)
            problem.AA{i,j} = spdiags(solin2{i}(:,j),0,size(solin2{i},1),size(solin2{i},1));
         end
         N_NT{i} =  speye(size(solin1{i},1));
         problem.BB{i} = solin1{i};
       end
       
       [sol2, err_reduc] = PGD_v6(problem,N_NT,options);
       for i=1:size(sol2.FF{1},2)
           sol2.FF{1}(:,i) = sol2.FF{1}(:,i)*sol2.alpha(i);
       end
       field = sol2.FF;
       display(['error reduction  : ' num2str(err_reduc(:)') ]);
       
    end
end

end

function Y = quotientEasy(U,V,options)
    %computes U/V in separated form
    
    %# of separated dimensions
    dim = numel(U);
    %# of dofs per dimension
    sz = cellfun(@(arg) size(arg,1),U);
    
    %# of terms in U & V
    NU = size(U{1},2);
    NV = size(V{1},2);
    
    Y = cell(size(U));
    R = cell(size(U));
    VY = cell(size(U));
    NVY = 0;

    ratio = 1;
    tol_enrich = options.res_reduc;
    tol_fp = options.fp_tol;
    max_enrichments = options.max_added_modes;
    max_fp = options.fp_max_iter;
    alphas = zeros(1,max_enrichments);
    %enrichment loop
    Nenrichments = 0;
    while (ratio>tol_enrich &&  Nenrichments<max_enrichments)
       fprintf(1,['E ' num2str(Nenrichments) ' ']);
       
        %initialization
        for d=1:dim
            R{d} = randn(sz(d),1);
            R{d} = R{d}/norm(R{d});
        end
        R_old = R;
        
        RtU = zeros(dim,NU);
        RtVY = zeros(dim,NVY);
        R2tV = zeros(dim,NV);
        
        for d=1:dim
            RtU(d,:) = R{d}'*U{d};
            if (Nenrichments>0)
                RtVY(d,:) = R{d}'*VY{d};
            end
            R2tV(d,:) = (R{d}').^2 * V{d};
        end
        
        diff_R = Inf;
        %fixed point loop
        Nfp = 0;
        while (Nfp<max_fp)
            fprintf(1,'.');
            for d = 1:dim
                mask = true(1,dim);
                mask(d) = false;
                %LHS:
                LHS = V{d}*prod(R2tV(mask,:),1)';
                
                %RHS: U contribution
                RHS = U{d}*prod(RtU(mask,:),1)';
                
                if Nenrichments>0
                %RHS: Y contribution
                RHS = RHS - VY{d}*prod(RtVY(mask,:),1)';
                end
                
                %solution
                R{d} = RHS./LHS;
                %normalization
                alpha_R = norm(R{d});
                R{d} = R{d}/alpha_R;
                %update of scalar quantities
                RtU(d,:) = R{d}'*U{d};
                if Nenrichments>0
                    RtVY(d,:) = R{d}'*VY{d};
                end
                R2tV(d,:) = (R{d}').^2 * V{d};
            end
            accum = 1;
            for d=1:dim
                accum = accum*R_old{d}'*R{d};
            end
            diff_R = sqrt(2-2*accum);
            R_old = R;
            Nfp = Nfp+1;
            
            if (diff_R < tol_fp) 
              fprintf(1,['CR ' num2str(ratio) ' \n']);
              break;
            end
        end
        if (Nfp==max_fp)
           fprintf(1,['CNR ' num2str(ratio) ' \n']);
        end
        
        R{1} = R{1}*alpha_R;
        %update Y and RHS
        for d = 1:dim
            Y{d} = [Y{d} R{d}];
            VY{d} = [VY{d} bsxfun(@times,V{d},R{d})];
        end
        NVY = NVY+NV;
        Nenrichments = Nenrichments+1;
        alphas(Nenrichments) = alpha_R;
        if Nenrichments==1
            alpha_1 = alpha_R;
        end
        ratio = alpha_R/alpha_1;
        

    end
    fprintf(1,'weight are : ');
    disp(num2str(alphas(1:Nenrichments)))
end


