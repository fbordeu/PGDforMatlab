%CREATE_FE_MESH
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Adrien Leygue (Adrien.Leygue@ec-nantes.fr)
%
function [X,INTERP,GEO] = CREATE_FE_MESH(mesh_type,varargin)
assert(ischar(mesh_type),'First argument should be a string');
mesh_type0= mesh_type;
mesh_type = upper(mesh_type);

nvar = numel(varargin);
switch mesh_type
    case 'GMSH'
        assert(nvar>=1,'GMSH requires a second argument');
        assert(ischar(varargin{1}),'second argument has to be the gmsh filename');
        [X,TOPO,NUMBER_TAG,PHYSICAL_TAG,~,PHYSICAL_NAMES,~,elem_names] = read_gmsh(varargin{1});
        assert(all(1:size(PHYSICAL_NAMES,1)==sort([PHYSICAL_NAMES{:,2}])),'Non continuous PHYSICAL TAGS');
        [sorted_tags,p] = sort([PHYSICAL_NAMES{:,2}]);
        PHYSICAL_NAMES = PHYSICAL_NAMES(p,3);
        
        if isempty(PHYSICAL_NAMES)
            tmp = [];
            for i=1:numel(PHYSICAL_TAG)
                tmp = [tmp;unique(PHYSICAL_TAG{i})]; %#ok<AGROW>
            end
            tmp = unique(tmp);
            for i=1:numel(tmp)
                PHYSICAL_NAMES{i} = ['Domain' num2str(i)];
            end
            sorted_tags = tmp;
        end
        for dd=1:numel(PHYSICAL_TAG)
            for tt = 1:numel(sorted_tags)
                PHYSICAL_TAG{dd}(PHYSICAL_TAG{dd}==sorted_tags(tt))=tt;
            end
        end
        
    case 'SEGMENT'
        assert(nvar==1,'SEGMENT requires a second argument');
        X = varargin{1};
        %validateattributes(X,{'double'},{'vector','finite'},'CREATE_FE_MESH','x',2);
        assert(numel(X)>1);
        X = X(:);
        Nx = numel(X);
        TOPO = {[1;Nx],[1:(Nx-1);2:Nx]'};
        NUMBER_TAG = {[1;2] (3:Nx+1)'};
        PHYSICAL_TAG = {[1;2] 3*ones(Nx-1,1)};
        PHYSICAL_NAMES = {'begin' 'end' 'Domain'};
        elem_names = {'POINT' 'EDGE_2'};
    case '2DPROD'
        assert(nvar~=3,'2DPROD requires exactly 3 argument ');
        x = varargin{1};
        y = varargin{2};
        Nx = numel(x);
        Ny = numel(y);
        [XX,YY] = meshgrid(x,y);
        X = [XX(:) YY(:)];
        id = reshape(1:(Nx*Ny),Ny,Nx);
        TOPO{1} = [id(1,1);id(1,end);id(end,end);id(end,1)];
        TOPO{2} = [id(1,1:end-1)' id(1,2:end)';...
            id(1:end-1,end) id(2:end,end);...
            id(end,end:-1:2)' id(end,end-1:-1:1)';...
            id(end:-1:2,1) id(end-1:-1:1,1)];
        TOPO{3} = zeros((Nx-1)*(Ny-1),4);
        TOPO{3}(:,1) = reshape(id(1:end-1,1:end-1),[],1);
        TOPO{3}(:,2) = reshape(id(1:end-1,2:end),[],1);
        TOPO{3}(:,3) = reshape(id(2:end,2:end),[],1);
        TOPO{3}(:,4) = reshape(id(2:end,1:end-1),[],1);
        
        
        NUMBER_TAG = {(1:4)'  4+(1:(2*(Nx-1)+2*(Ny-1)))' (4+2*(Nx-1)+2*(Ny-1))+(1:(Nx-1)*(Ny-1))'};
        PHYSICAL_TAG = {[1;2;3;4] [ones(Nx-1,1)*5; ones(Ny-1,1)*6;ones(Nx-1,1)*7;ones(Ny-1,1)*8]   [9*ones((Ny-1)*(Nx-1),1)]};
        PHYSICAL_NAMES = {'bl' 'br' 'ul' 'ur' 'down' 'right' 'up' 'left' 'Domain'};
        elem_names = {'POINT' 'EDGE_2' 'QUAD_4'};
    case '2DPRODDISCRETE'
        assert(nvar~=3,'2DPROD requires exactly 3 argument ');
        x = varargin{1};
        y = varargin{2};
        Nx = numel(x);
        Ny = numel(y);
        [XX,YY] = meshgrid(x,y);
        X = [XX(:) YY(:)];
        
        TNN = Nx*Ny;% Total Number Of Nodes;
        
        id = reshape(1:TNN,Ny,Nx);
        
        TOPO{1} = [id(:)];
        TOPO{2} = [id(1,1:end-1)' id(1,2:end)';...
            id(1:end-1,end) id(2:end,end);...
            id(end,end:-1:2)' id(end,end-1:-1:1)';...
            id(end:-1:2,1) id(end-1:-1:1,1)];
        TOPO{3} = zeros((Nx-1)*(Ny-1),4);
        TOPO{3}(:,1) = reshape(id(1:end-1,1:end-1),[],1);
        TOPO{3}(:,2) = reshape(id(1:end-1,2:end),[],1);
        TOPO{3}(:,3) = reshape(id(2:end,2:end),[],1);
        TOPO{3}(:,4) = reshape(id(2:end,1:end-1),[],1);
        
        
        NUMBER_TAG = {(1:TNN)'   TNN+(1:(2*(Nx-1)+2*(Ny-1)))' (TNN+2*(Nx-1)+2*(Ny-1))+(1:(Nx-1)*(Ny-1))'};
        PHYSICAL_TAG = {[ones(TNN,1)] [ones(Nx-1,1)*2; ones(Ny-1,1)*3;ones(Nx-1,1)*4;ones(Ny-1,1)*5]   [6*ones((Ny-1)*(Nx-1),1)]};
        PHYSICAL_NAMES = {'NODAL' 'down' 'right' 'up' 'left' 'Domain'};
        elem_names = {'POINT' 'EDGE_2' 'QUAD_4'};
    case '3DPROD'
        assert(nvar~=4,'3DPROD requires exactly 4 argument ');
        x = varargin{1};
        y = varargin{2};
        z = varargin{3};
        Nx = numel(x);
        Ny = numel(y);
        Nz = numel(z);
        [XX,YY,ZZ] = meshgrid(x,y,z);
        X = [XX(:) YY(:) ZZ(:)];
        id = reshape(1:(Nx*Ny*Nz),Ny,Nx,Nz);
        %% Points
        TOPO{1} = [id(1,1,1);id(1,end,1);id(end,end,1);id(end,1,1);...
            id(1,1,end);id(1,end,end);id(end,end,end);id(end,1,end)];
        %% edges
        TOPO{2} = [id(1,1:end-1,1)' id(1,2:end,1)';...% bottom
            id(1:end-1,end,1) id(2:end,end,1);...
            id(end,end:-1:2,1)' id(end,end-1:-1:1,1)';...
            id(end:-1:2,1,1) id(end-1:-1:1,1,1);...
            id(1,1:end-1,end)' id(1,2:end,end)';...% top
            id(1:end-1,end,end) id(2:end,end,end);...
            id(end,end:-1:2,end)' id(end,end-1:-1:1,end)';...
            id(end:-1:2,1,end) id(end-1:-1:1,1,end);...
            %%TODO the Vertical EDGEs
            squeeze(id(1,1,1:end-1)) squeeze(id(1,1,2:end));...% verticals
            squeeze(id(end,1,1:end-1)) squeeze(id(end,1,2:end));...
            squeeze(id(end,end,1:end-1)) squeeze(id(end,end,2:end));...
            squeeze(id(1,end,1:end-1)) squeeze(id(1,end,2:end));....
            ];
        
        %%TODO to code the ceil and the vertical faces
        %             TOPO{3} = zeros((Nx-1)*(Ny-1),4);
        %             TOPO{3}(:,1) = reshape(id(1:end-1,1:end-1,1),[],1);
        %             TOPO{3}(:,2) = reshape(id(1:end-1,2:end,1),[],1);
        %             TOPO{3}(:,3) = reshape(id(2:end,2:end,1),[],1);
        %             TOPO{3}(:,4) = reshape(id(2:end,1:end-1,1),[],1);
        TOPO{3} = zeros(2*(Nx-1)*(Ny-1)+2*(Nx-1)*(Nz-1)+2*(Ny-1)*(Nz-1),4);
        TOPO{3}(:,1) = [reshape(id(1:end-1,1:end-1,1),[],1); reshape(id(1:end-1,1:end-1,end),[],1);...
            reshape(id(1,1:end-1,1:end-1),[],1); reshape(id(end,1:end-1,1:end-1),[],1);...
            reshape(id(1:end-1,1,1:end-1),[],1); reshape(id(1:end-1,end,1:end-1),[],1);...
            ];
        TOPO{3}(:,2) = [reshape(id(1:end-1,2:end,1),[],1); reshape(id(1:end-1,2:end,end),[],1);...
            reshape(id(1,2:end,1:end-1),[],1); reshape(id(end,2:end,1:end-1),[],1);...
            reshape(id(2:end,1,1:end-1),[],1); reshape(id(2:end,end,1:end-1),[],1);...
            ];
        TOPO{3}(:,3) = [reshape(id(2:end,2:end,1),[],1); reshape(id(2:end,2:end,end),[],1);...
            reshape(id(1,2:end,2:end),[],1); reshape(id(end,2:end,2:end),[],1);...
            reshape(id(2:end,1,2:end),[],1); reshape(id(2:end,end,2:end),[],1);...
            ];
        TOPO{3}(:,4) = [reshape(id(2:end,1:end-1,1),[],1); reshape(id(2:end,1:end-1,end),[],1);...
            reshape(id(1,1:end-1,2:end),[],1); reshape(id(end,1:end-1,2:end),[],1);...
            reshape(id(1:end-1,1,2:end),[],1); reshape(id(1:end-1,end,2:end),[],1);...
            ];
        
        TOPO{4} = zeros((Nx-1)*(Ny-1)*(Nz-1),8);
        TOPO{4}(:,1) = reshape(id(1:end-1,1:end-1,1:end-1),[],1);
        TOPO{4}(:,2) = reshape(id(1:end-1,2:end  ,1:end-1),[],1);
        TOPO{4}(:,3) = reshape(id(2:end  ,2:end  ,1:end-1),[],1);
        TOPO{4}(:,4) = reshape(id(2:end  ,1:end-1,1:end-1),[],1);
        TOPO{4}(:,5) = reshape(id(1:end-1,1:end-1,2:end  ),[],1);
        TOPO{4}(:,6) = reshape(id(1:end-1,2:end  ,2:end  ),[],1);
        TOPO{4}(:,7) = reshape(id(2:end  ,2:end  ,2:end  ),[],1);
        TOPO{4}(:,8) = reshape(id(2:end  ,1:end-1,2:end  ),[],1);
        
        %             NUMBER_TAG = {(1:8)'  8+(1:(4*(Nx-1)+4*(Ny-1)))' (8+4*(Nx-1)+4*(Ny-1))+(1:(Nx-1)*(Ny-1))' (8+4*(Nx-1)+4*(Ny-1)+(Nx-1)*(Ny-1))+(1:(Nx-1)*(Ny-1)*(Nz-1))' };
        NUMBER_TAG = {(1:8)'  8+(1:(4*(Nx-1)+4*(Ny-1)+4*(Nz-1)))' (8+4*(Nx-1)+4*(Ny-1)+4*(Nz-1))+(1:(2*((Nx-1)*(Ny-1))+2*((Nx-1)*(Nz-1))+2*((Ny-1)*(Nz-1))))' (8+4*(Nx-1)+4*(Ny-1)+4*(Nz-1)+2*(Nx-1)*(Ny-1)+2*(Nx-1)*(Nz-1)+2*(Ny-1)*(Nz-1))+(1:(Nx-1)*(Ny-1)*(Nz-1))' };
        PHYSICAL_TAG = {[1;2;3;4;5;6;7;8]...
            [ones(Nx-1,1)*9; ones(Ny-1,1)*10;ones(Nx-1,1)*11;ones(Ny-1,1)*12;ones(Nx-1,1)*13; ones(Ny-1,1)*14;ones(Nx-1,1)*15;ones(Ny-1,1)*16;ones(Nz-1,1)*17;ones(Nz-1,1)*18;ones(Nz-1,1)*19;ones(Nz-1,1)*20;]...
            [21*ones((Ny-1)*(Nx-1),1); 22*ones((Ny-1)*(Nx-1),1); 23*ones((Nx-1)*(Nz-1),1); 24*ones((Nx-1)*(Nz-1),1); 25*ones((Ny-1)*(Nz-1),1); 26*ones((Ny-1)*(Nz-1),1);]...
            [27*ones((Ny-1)*(Nx-1)*(Nz-1),1)]};
        PHYSICAL_NAMES = {'bbb' 'btb' 'ttb' 'tbb' 'bbt' 'btt' 'ttt' 'tbt'  'bdown' 'bright' 'bup' 'bleft'   'tdown' 'tright' 'tup' 'tleft'   'bb_bt' 'tb_bt' 'tt_bt' 'bt_bt'   'floor' 'ceiling' 'left' 'right' 'back' 'front' 'Domain'};
        flip_dom = find(ismember(PHYSICAL_NAMES,{'floor','right','back'}));
        mask = ismember(PHYSICAL_TAG{3},flip_dom);
        TOPO{3}(mask,:) = fliplr(TOPO{3}(mask,:));
        elem_names = {'POINT' 'EDGE_2' 'QUAD_4' 'HEX_8'};
    case '3DPRODDISCRETE'
        assert(nvar~=4,'3DPROD requires exactly 4 argument ');
        x = varargin{1};
        y = varargin{2};
        z = varargin{3};
        Nx = numel(x);
        Ny = numel(y);
        Nz = numel(z);
        [XX,YY,ZZ] = meshgrid(x,y,z);
        X = [XX(:) YY(:) ZZ(:)];
        id = reshape(1:(Nx*Ny*Nz),Ny,Nx,Nz);
        %% Points
        TOPO{1} = [ id(:) ];
        %% edges
        TOPO{2} = [id(1,1:end-1,1)' id(1,2:end,1)';...% bottom
            id(1:end-1,end,1) id(2:end,end,1);...
            id(end,end:-1:2,1)' id(end,end-1:-1:1,1)';...
            id(end:-1:2,1,1) id(end-1:-1:1,1,1);...
            id(1,1:end-1,end)' id(1,2:end,end)';...% top
            id(1:end-1,end,end) id(2:end,end,end);...
            id(end,end:-1:2,end)' id(end,end-1:-1:1,end)';...
            id(end:-1:2,1,end) id(end-1:-1:1,1,end)];....
            %%TODO the Vertical EDGEs
        %id(1,1,1:end-1)' id(1,1,2:end)';...% verticals
        %id(end,1,1:end-1) id(end,1,2:end);...
        %id(end,end:-1:2,end)' id(end,end-1:-1:1,end)';...
        %id(end:-1:2,1,end) id(end-1:-1:1,1,end);....
        %];
        
        %%TODO to code the ceil and the vertical faces
        TOPO{3} = zeros((Nx-1)*(Ny-1),4);
        TOPO{3}(:,1) = reshape(id(1:end-1,1:end-1,1),[],1);
        TOPO{3}(:,2) = reshape(id(1:end-1,2:end,1),[],1);
        TOPO{3}(:,3) = reshape(id(2:end,2:end,1),[],1);
        TOPO{3}(:,4) = reshape(id(2:end,1:end-1,1),[],1);
        
        TOPO{4} = zeros((Nx-1)*(Ny-1)*(Nz-1),8);
        TOPO{4}(:,1) = reshape(id(1:end-1,1:end-1,1:end-1),[],1);
        TOPO{4}(:,2) = reshape(id(1:end-1,2:end  ,1:end-1),[],1);
        TOPO{4}(:,3) = reshape(id(2:end  ,2:end  ,1:end-1),[],1);
        TOPO{4}(:,4) = reshape(id(2:end  ,1:end-1,1:end-1),[],1);
        TOPO{4}(:,5) = reshape(id(1:end-1,1:end-1,2:end  ),[],1);
        TOPO{4}(:,6) = reshape(id(1:end-1,2:end  ,2:end  ),[],1);
        TOPO{4}(:,7) = reshape(id(2:end  ,2:end  ,2:end  ),[],1);
        TOPO{4}(:,8) = reshape(id(2:end  ,1:end-1,2:end  ),[],1);
        
        TNN = Nx*Ny*Nz;% Total Number Of Nodes;
        
        NUMBER_TAG = {(1:TNN)'  TNN+(1:(4*(Nx-1)+4*(Ny-1)))' (TNN+4*(Nx-1)+4*(Ny-1))+(1:(Nx-1)*(Ny-1))' (TNN+4*(Nx-1)+4*(Ny-1)+(Nx-1)*(Ny-1))+(1:(Nx-1)*(Ny-1)*(Nz-1))' };
        PHYSICAL_TAG = {[ones(TNN,1)]...
            [ones(Nx-1,1)*2; ones(Ny-1,1)*3;ones(Nx-1,1)*4;ones(Ny-1,1)*5;ones(Nx-1,1)*6; ones(Ny-1,1)*7;ones(Nx-1,1)*8;ones(Ny-1,1)*9;]...
            [10*ones((Ny-1)*(Nx-1),1) ]...
            [11*ones((Ny-1)*(Nx-1)*(Nz-1),1)]};
        PHYSICAL_NAMES = {'NODAL'  'bdown' 'bright' 'bup' 'bleft'   'tdown' 'tright' 'tup' 'tleft'  'floor'  'Domain'};
        elem_names = {'POINT' 'EDGE_2' 'QUAD_4' 'HEX_8'};
    case 'DISCRETE'
        %             assert(nvar==1,'DISCRETE requires a second argument');
        %             X = varargin{1};
        %             validateattributes(X,{'double'},{'vector','finite'},'CREATE_FE_MESH','x',2);
        %             Nx = numel(X);
        %             TOPO = {(1:Nx)'};
        %             NUMBER_TAG = {(1:Nx)'};
        %             PHYSICAL_TAG = {ones(Nx,1)};
        %             PHYSICAL_NAMES = {'Discrete'};
        %             elem_names = {'EDGE_1'};
        assert(nvar==1,'DISCRETE requires a second argument');
        X = varargin{1};
        validateattributes(X,{'double'},{'vector','finite'},'CREATE_FE_MESH','x',2);
        assert(numel(X)>1);
        X = X(:);
        Nx = numel(X);
        TOPO = {(1:Nx)',[1:(Nx-1);2:Nx]'};
        NUMBER_TAG = {(1:Nx)' Nx+(1:(Nx-1))'};
        PHYSICAL_TAG = {ones(Nx,1) 2*ones(Nx-1,1)};
        PHYSICAL_NAMES = {'Discrete' 'Continuous'};
        elem_names = {'POINT' 'EDGE_2'};
    otherwise
        error(['Don''t know how to treat : ' mesh_type0]);
end

GEO = GEOMETRY(TOPO,NUMBER_TAG,PHYSICAL_TAG,PHYSICAL_NAMES,elem_names);
options = {};
for i=1:numel(PHYSICAL_NAMES)
    options = [options { 'SUB' PHYSICAL_NAMES{i} }]; %#ok<AGROW>
end
INTERP = FEM_INTERP(GEO,'ALL','DEFAULT',options);

end