%FEM_FIELD CLASS
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Adrien Leygue (Adrien.Leygue@ec-nantes.fr)
%
classdef FEM_FIELD < handle
    properties
        name = '';
        Ncomp = 1;
        Nmodes = 1;
        Nnodes = 0;
        N_PHI = 0;
        CURRENT_INTERP = [];
    %end
    %properties (Access='private', Hidden=true)
        type = 0;
        INTERP = [];
        VALUES = [];
        PARENT_FIELD = [];
        FCT = [];
        CURRENT_ID = 0;
        
        OFFSET_LOC = 0;
        OFFSET_GLOB = 0;
        N_LOC = 0;
        N_GLOB = 0;
        
    end
    methods
        function obj = FEM_FIELD(name,a,b,Ncomp,Nmodes)
            assert(ischar(name),'First argument should be a string');
            if nargin==0
                return;
            end
            obj.name = name;
            
            if isa(a,'FEM_INTERP')
                obj.type=1;
                obj.INTERP = a;
                obj.Nnodes = obj.INTERP(1).MAX_VALUES;
                switch nargin
                    case 2
                        b = [];
                        Ncomp = 1;
                        Nmodes = 1;
                    case 3
                        Ncomp = size(b,2);
                        Nmodes = 1;
                    case 4
                        Nmodes = 1;
                        if ~isempty(b)
                            Nmodes = round(size(b,2)/Ncomp);
                            assert((Nmodes*Ncomp)==size(b,2),'size(Values,2) should be a multiple of the number of components');
                        end
                    case 5
                        if ~isempty(b)
                            assert((Nmodes*Ncomp)==size(b,2),'size(Values,2) should be a Nmodes*Ncomp');
                        end
                end
                
                obj.VALUES = b;
                obj.Ncomp = Ncomp;
                obj.Nmodes = Nmodes;
            elseif isa(a,'FEM_FIELD')
                obj.type = 2;
                obj.PARENT_FIELD = a;
                obj.FCT = b;
                switch nargin
                    case 3
                        Ncomp = 1;
                        Nmodes = 1;
                    case 4
                        Nmodes = 1;
                end
                obj.Ncomp = Ncomp;
                obj.Nmodes = Nmodes;
            else
                error('Invalid first argument');
            end
        end
        
        function result = eval_at_integration_points(obj,domain,integ_param)
             interp_id = obj.INTERP.find_domain(domain);
             assert(interp_id>0,'Unknown domain');
             pc_id = (integ_param==obj.INTERP(interp_id).PCIntegrationParam);
             assert(nnz(pc_id)==1,'no precomputed gauss values');
             assert(~isempty(obj.VALUES),'Field has no nodal values');
             result = obj.INTERP(interp_id).PCMatrix{pc_id}*obj.VALUES(:,1:obj.Ncomp);    
        end
        
        function result = get_dofs(obj,varargin)
            %gives the numbering of the d.o.f.s of a field in a vector of fields
            %Valid syntax
            %
            %OBJ.get_dofs(FIELD_NAME,DOMAIN_NAME,COMPONENTS);
            % FIELDS.get_dofs('U','BND',[1 3]) = dof number of the
            % components 1 and 3 of the field 'U' on the domain 'BND'
            %
            %OBJ.get_dofs(DOMAIN_NAME,COMPONENTS);
            %FIELD_NAME may be ommited if numel(OBJ)==1
            %
            %OBJ.get_dofs(FIELD_NAME, DOMAIN_NAME);
            %COMPONENTS may be ommited. Default value= 1:Number_of_components
            %
            %OBJ.get_dofs(DOMAIN_NAME);
            %COMPONENTS *AND* FIELD may be ommited.
            
            %process the inputs
            [field_id,dom_id,~,comp_id] = obj.get_ids(varargin);
            
            %identifies the dofs numbers from the connectivity of the
            %domain and the components
            result = bsxfun(@plus,obj(field_id).INTERP(dom_id).get_dofs(), (comp_id-1)*obj(field_id).Nnodes);
            %compute the offset if numel(obj) >1
            offset = 0;
            for i=1:(field_id-1)
                offset = offset + obj(i).Nnodes*obj(i).Ncomp;
            end
            result = result(:) + offset;
        end
        
        function result = get_values(obj,varargin)
            %gives the values of the d.o.f.s of a field in a vector of fields
            %Valid syntax
            %
            %U = OBJ.get_values(FIELD_NAME,DOMAIN_NAME,COMPONENTS);
            %U = FIELDS.get_values('U','BND',[1 3]) = dof values of the
            % components 1 and 3 of the field 'U' on the domain 'BND'
            %
            %U = OBJ.get_values(DOMAIN_NAME,COMPONENTS);
            %FIELD_NAME may be ommited if numel(OBJ)==1
            %
            %U = OBJ.get_values(FIELD_NAME, DOMAIN_NAME);
            %COMPONENTS may be ommited. Default value= 1:Number_of_components
            %
            %U = OBJ.get_values(DOMAIN_NAME);
            %COMPONENTS *AND* FIELD may be ommited.
            %
            %NB: It is possible to store several values simultaneously.
            %each value is called a mode.
            %If two numeric inputs are provided, the second is the mode
            %numbers: 
            %U = OBJ.get_values(FIELD_NAME,DOMAIN_NAME,COMPONENTS,MODES);
            %U = FIELDS.get_values('U','BND',[1 3],[4 5]) = dof values of the
            % components 1 and 3 of modes 4 and 5 the field 'U' on the domain 'BND'
            
            [field_id,dom_id,cols] = obj.get_ids(varargin);
            result = obj(field_id).VALUES(obj(field_id).INTERP(dom_id).get_dofs(),cols);
        end
        
        function set_values(obj,values,varargin)
            %Sets the values of the d.o.f.s of a field in a vector of fields
            %Valid syntax
            %
            %OBJ.set_values(Values,FIELD_NAME,DOMAIN_NAME,COMPONENTS);
            %U = FIELDS.get_values(Values,'U','BND',[1 3])sets the values of the
            % components 1 and 3 of the field 'U' on the domain 'BND' to
            % values
            %
            %OBJ.get_values(Values,DOMAIN_NAME,COMPONENTS);
            %FIELD_NAME may be ommited if numel(OBJ)==1
            %
            %OBJ.set_values(Values,FIELD_NAME, DOMAIN_NAME);
            %COMPONENTS may be ommited. Default value= 1:Number_of_components
            %
            %OBJ.set_values(Values,DOMAIN_NAME);
            %COMPONENTS *AND* FIELD may be ommited.
            %
            %NB: It is possible to store several values simultaneously.
            %each value is called a mode.
            %If two numeric inputs are provided, the second is the mode
            %numbers: 
            %OBJ.set_values(Values,FIELD_NAME,DOMAIN_NAME,COMPONENTS,MODES);
            %FIELDS.set_values(Values,'U','BND',[1 3],[4 5])  sets the values of the
            % components 1 and 3 of modes 4 and 5 the field 'U' on the domain 'BND'
            
            [field_id,dom_id,cols,~,mode_id] = obj.get_ids(varargin);
            lines = obj(field_id).INTERP(dom_id).get_dofs;
            obj(field_id).Nmodes = max([obj(field_id).Nmodes max(mode_id)]);
            obj(field_id).VALUES(lines,cols) = values; %#ok<NASGU>
        end
        
        function set_all_values(obj,values)
            
            if iscell(values)
                for i= 1:numel(values)
                    obj(i).set_all_values(values{i});
                end
                return
            end
                
            offsets  = [obj.Ncomp].*[obj.Nnodes];
            offsets = [0;cumsum(offsets(:))];
            offsets=offsets(1:(end-1));
            for i=1:numel(obj)    
                if (~isscalar(values))
                obj(i).Nmodes = size(values,2);
                tmp = values(offsets(i)+(1:(obj(i).Nnodes*obj(i).Ncomp)),:);
                else
                   tmp = values*ones(obj(i).Nnodes,obj(i).Ncomp*obj(i).Nmodes);
                end
                obj(i).VALUES = reshape(tmp,obj(i).Nnodes,obj(i).Ncomp*obj(i).Nmodes);
            end
        end
        
        function result = get_all_values(obj)
            result = reshape(obj(1).VALUES,obj(1).Nnodes*obj(1).Ncomp,obj(1).Nmodes);
            for i=2:numel(obj)
                result = [result;reshape(obj(i).VALUES,obj(i).Nnodes*obj(i).Ncomp,obj(i).Nmodes)];
            end
        end
        
        
        function set_modes_number(obj,N)
            if (obj.Nmodes ==N), return; end
            tmp = obj.Nmodes;
            
            obj.Nmodes = N;
            if isempty(obj.VALUES)
                return;
            end
            if N<tmp
                obj.VALUES = obj.VAUES(:,1:(obj.Ncomp*obj.Nmodes));
            else
                obj.VALUES(end,obj.Ncomp*obj.Nmodes) = 0;
            end
            
        end
        
        function Nelem = set_current_domain(obj,name)
            %Nelem = OBJ.set_current_domain(DOMAIN)
            %For all fields in OBJ, identifies if the field has a specific mesh
            %for the domain DOMAIN.
            %If yes, sets a shortcut to access this domain faster.
            %Nelem is the number of elements on the domain
            
            Nelem = 0;
            %loop over the fields
            for i=1:numel(obj)
                %different treatment id the field is a FE field or a
                %function of a FE field
                switch obj(i).type
                    case 1
                        obj(i).CURRENT_ID = obj(i).INTERP.find_domain(name);
                        if obj(i).CURRENT_ID
                            obj(i).CURRENT_INTERP = obj(i).INTERP(obj(i).CURRENT_ID);
                            Nelem = obj(i).CURRENT_INTERP.N_ELEMENTS;                       
                        else
                            obj(i).CURRENT_INTERP = [];
                        end
                    case 2
                            obj(i).PARENT_FIELD.set_current_domain(name);
                end
            end
        end
        
        function [Nglob,Nloc] = compute_offsets(obj)
            %[Nglob,Nloc] = OBJ.compute_offsets();
            %compute local and global offsets for the dofs of all fields of OBJ
            %
            %Global offset is the offset between the numbering in the
            %connectivity table and the numbering in the assembled system
            %
            %Local offset is the offset of the dofs position in the local
            %matrix (during assembly)
            
            assert(all([obj.type]==1),'invalid FEM_FIELD for compute_offsets');
            Nloc = 0;
            Nglob = 0;
            for i=1:numel(obj)
                obj(i).OFFSET_LOC = Nloc;
                obj(i).OFFSET_GLOB = Nglob;
                Nglob = Nglob + obj(i).Nnodes*obj(i).Ncomp;

                if obj(i).CURRENT_ID
                    obj(i).N_PHI = obj(i).CURRENT_INTERP.ELEM.Nphi;
                else
                    obj(i).N_PHI =0;
                end
                
                Nloc = Nloc+ obj(i).N_PHI*obj(i).Ncomp;
                
            end
            for i=1:numel(obj)
                obj(i).N_LOC = Nloc;
                obj(i).N_GLOB = Nglob;
            end
        end
        
        
        function result = get_current_elements(obj)
            %result = OBJ.get_current_elements()
            %returns the vector of elements of all the fields for the current 
            %domain (defined through OBJ.set_current_domain)
            %if a field has no mesh for the current domain, an empty
            %element structure is given.
            result = ELEMENT;
            result = repmat(result,size(obj));
            for i=1:numel(obj)
                if obj(i).type == 1
                    if obj(i).CURRENT_ID
                    result(i) = obj(i).CURRENT_INTERP.ELEM;
                    end
                else
                    result(i) = obj(i).PARENT_FIELD.get_current_elements();
                end
            end
        end
        
        function result = get_local_index(obj)
            %This function computes the assembly position in the local
            %matrix of all the components of the field OBJ
            %NB: obj cannot be a vector of fields
            result = zeros(obj.Ncomp,obj.N_PHI);
            for i=1:obj.Ncomp
                result(i,:) = (1:obj.N_PHI)+obj.OFFSET_LOC+(i-1)*obj.N_PHI;
            end
        end
        
        function result = get_detailled_map(obj,el,comp)
            %This function computes the detailled assembly position, of a
            %for the shape functions of field obj, component comp in
            %element el.
            result = obj.CURRENT_INTERP.get_element_dofs(el)+obj.CURRENT_INTERP.MAX_VALUES*(comp-1)+ obj.OFFSET_GLOB;
        end
        
        function result = get_local_to_global_map(obj,el)
            %This function computes the local to global mapping for element el for all
            %fields in obj.
%             result = zeros(1,obj(1).N_LOC);
%             for i=1:numel(obj)
%                 if (obj(i).CURRENT_ID)
%                 for j=1:obj(i).Ncomp
%                     result(obj(i).OFFSET_LOC+(j-1)*obj(i).N_PHI + (1:obj(i).N_PHI)) = obj(i).CURRENT_INTERP.get_element_dofs(el)+obj(i).CURRENT_INTERP.MAX_VALUES*(j-1)+ obj(i).OFFSET_GLOB;
%                 end
%                 end
%             end
            result = [];
%             for i=1:numel(obj)
%                 if (obj(i).CURRENT_ID)
%                 for j=1:obj(i).Ncomp
%                     result = [ result obj(i).CURRENT_INTERP.get_element_dofs(el)+(obj(i).CURRENT_INTERP.MAX_VALUES*(j-1)+ obj(i).OFFSET_GLOB)];
%                 end
%                 end
%             end
             for i=1:numel(obj)
                 if (obj(i).CURRENT_ID)
                        tmp = bsxfun(@plus,obj(i).OFFSET_GLOB+obj(i).CURRENT_INTERP.get_element_dofs(el)',obj(i).Nnodes*(0:obj(i).Ncomp-1));
                        result = [result tmp(:)']; %#ok<AGROW>
                 end
             end

        end
        
        function result = get_field_connectivity(obj)
            assert(numel(obj)==1,'This function is only applicable to a single field');
            result = cell(1,obj.Ncomp);
            for c = 1:obj.Ncomp
            result{c} = obj.CURRENT_INTERP.CONNECTIVITY+obj.OFFSET_GLOB+(c-1)*obj.Nnodes;
            end
        end
        
        function result = eval_phi_integration(obj,el,J,ders,comp) %#ok<INUSL,INUSD>
%             dim = obj.CURRENT_INTERP.ELEM.dim;
%             if(el ~=obj.last_el)
%                 obj.buffer_phi = cell(obj.CURRENT_INTERP.ELEM.Neval,2^dim);
%                 obj.buffer_phi_mask = false(1,2^dim);
%                 obj.last_el = el;
%             end
%             %tmp = num2cell(ders+1);
%             idx = 1+ (2.^(0:(dim-1)))*(ders);
%             %idx2 = sub2ind(2*ones(1,dim),tmp{:});
%             if obj.buffer_phi_mask(idx)
%                 result = obj.buffer_phi{idx};
%             else
               result = obj.CURRENT_INTERP.ELEM.gen_eval_phi(J,ders);
%              obj.buffer_phi{idx} = result;
%              obj.buffer_phi_mask(idx) = 1;
%             end
        end
        function result = eval_field_integration(obj,el,J,ders,comp,mode)
            %evaluate the comp^th component of the mode^th mode of field
            %obj at its integration points in element el.
            %a derivative identified by ders might be taken. The jacobian
            %matrix id J (cell array of matrices, one for each integration point)
            if nargin<6
                mode=1;
            end
            if nargin<5
                comp=1:obj.Ncomp;
            end
            if isempty(comp)
                comp=1:obj.Ncomp;
            end
            switch obj.type
                case 1
                    Fi = obj.VALUES(obj.CURRENT_INTERP.get_element_dofs(el),comp+(mode-1)*(obj.Ncomp));
                    result = obj.CURRENT_INTERP.ELEM.eval_field(Fi,J,ders);
                    return;
                case 2
                    result = obj.PARENT_FIELD.eval_field_integration(el,J,[],[],mode);
                    result = obj.FCT(result');
                    result = result(:,comp)';
                    return;
            end
        end
        
        function W = prepare_integration(obj,varargin)
            %Prepare the elements of the current domain for integration
            %(evaluation of the shape functions and their gradients at the
            %integration point)
            %extra arguents are passed to the ELEMENT.
            for i=1:numel(obj)
                switch obj(i).type
                    case 1
                        if obj(i).CURRENT_ID
                            W = obj(i).CURRENT_INTERP.ELEM.prepare_integration(varargin{:});
                        end
                    case 2
                        W = obj(i).PARENT_FIELD.prepare_integration(varargin{:});
                end
            end
        end
        
        function [J,detJ,iJ] = jacobian_integration(obj,el)
            %Computes the jacobian matrix at the integration points
            Fi = obj.VALUES(obj.CURRENT_INTERP.get_element_dofs(el),:);
            [J,detJ] = obj.CURRENT_INTERP.ELEM.jacobian(Fi);
            if nargin>=3
                iJ = cellfun(@inv,J,'uniformoutput','false');
            end
        end
        
        function LSfit(obj,integ_param,integ_values)
            M = obj.INTERP(1).PCMatrix{find(integ_param==obj.INTERP(1).PCIntegrationParam),1};
            W = obj.INTERP(1).PCIntegrationWeights{find(integ_param==obj.INTERP(1).PCIntegrationParam)};
            NW = numel(W);
            W = spdiags(W(:),0,NW,NW);
            values = (M'*W*M)\(M'*W*integ_values);
           obj.set_all_values(values);
        end
        
        function precompute_interpolants(obj,domain,integ_param,all_fields)
            if isempty(domain)
                return;
            end
            if iscellstr(domain)
                for i=1:numel(domain)
                    precompute_interpolants(obj,domain{i},integ_param,all_fields)
                end
            end
            all_interp = obj.INTERP(:);
            for i=1:numel(all_fields)
                all_interp = [all_interp;all_fields(i).INTERP(:)] ;%#ok<AGROW>
            end
            all_interp = unique(all_interp);
            all_interp = all_interp(all_interp.find_domain(domain));
            Ninterp = numel(all_interp);
            map_interp = all_interp((ismember(all_interp,obj.INTERP)));
            loc_id = zeros(1,Ninterp);
            %attention: virer les interps d�l� pr�calcul�s pour ce
            %param�tre
            mask = zeros(1,Ninterp);
            for i=1:Ninterp
                mask(i) = any(all_interp(i).PCIntegrationParam==integ_param);
            end
            all_interp = all_interp(~mask);
            Ninterp = numel(all_interp);
            if Ninterp==0;
                return;
            end
            
            for i=1:Ninterp
                    all_interp(i).PCIntegrationParam = [all_interp(i).PCIntegrationParam integ_param];
                    loc_id(i) = numel(all_interp(i).PCIntegrationParam);
            end
            all_elements = vertcat(all_interp.ELEM);
            all_elements = unique(all_elements);
            for e = 1:numel(all_elements)
                elem_weights = all_elements(e).prepare_integration(integ_param);
            end
            elem_weights = elem_weights(:);
            map_element = map_interp.ELEM;
            Nelem = map_interp.N_ELEMENTS;
            dim = map_element.dim;
            Ninteg = numel(elem_weights);
            for i=1:Ninterp
                all_interp(i).NPCIntegrationPoints(loc_id(i)) = Nelem*Ninteg;
                all_interp(i).NPCIntegrationPointsPerElement(loc_id(i)) = Ninteg;
                
                all_interp(i).PCIntegrationWeights{loc_id(i)} = zeros(Nelem*Ninteg,1);
                for d = 1:(dim+1)
                    all_interp(i).PCMatrix{loc_id(i),d} = zeros(Nelem*Ninteg,all_interp(i).ELEM.Nphi);
                end
            end
            
            for e = 1:Nelem
                Fi = obj.VALUES(map_interp.CONNECTIVITY(e,:),:);
                [J,detJ] = map_element.jacobian(Fi);
                tmp = (1:Ninteg)+(e-1)*Ninteg;
                for i = 1:Ninterp
                    all_interp(i).PCIntegrationWeights{loc_id(i)}(tmp) = detJ.*elem_weights;
                    all_interp(i).PCMatrix{loc_id(i),1}(tmp,:) = all_interp(i).ELEM.eval_phi;
                    for p = 1:Ninteg
                        dN = J{p} \  all_interp(i).ELEM.eval_dphidxi{p};
                    for d = 1:dim
                        all_interp(i).PCMatrix{loc_id(i),d+1}(p+(e-1)*Ninteg,:) = dN(d,:)';
                    end
                    end
                end
            end
            for i=1:Ninterp
                tmp = kron(all_interp(i).CONNECTIVITY,ones(Ninteg,1));
                tmp2 = repmat((1:Nelem*Ninteg)',[1 all_interp(i).ELEM.Nphi]);
                for d = 1:(dim+1)
                    all_interp(i).PCTableValues{loc_id(i),d} = all_interp(i).PCMatrix{loc_id(i),d};
                    all_interp(i).PCTableIndex{loc_id(i),d} = tmp;                    
                    all_interp(i).PCMatrix{loc_id(i),d} = sparse(tmp2(:),tmp(:),all_interp(i).PCMatrix{loc_id(i),d}(:),Nelem*Ninteg,all_interp(i).MAX_VALUES);
                end
            end
            
        end
        
        function data = export_xdmf(obj,filename_data,name,varargin)
            %XDMF /PXDMF export function
            %Usage
            %
            %COORDS.export_xdmf(FILENAME,DOMAIN,FIELD1,FIELD2,...)
            % export the domain DOMAIN described by the coordinate Field COORDS to xdmf.
            % Additionally exports fields FIELD1, FIELD2,... on that domain.
            %
            %data = COORDS.export_xdmf('',DOMAIN,FIELD1,FIELD2,...)
            % same but the file is not written. It can be later written
            % with writepxdmf(data) (don't forget to set data.filename!!)
            %
            %data = COORDS_X.export_xdmf('',DOMAIN_D1,FIELD1,FIELD2,...)
            %data = COORDS_Y.export_xdmf(data,DOMAIN_D2,FIELD1,FIELD2,...)
            %writepxdmf(data)
            %Writes the pxdmf file of the separated COORDS_X * COORDS_Y
            %representation of the different fields.
            
            
            assert(obj.type==1,'THIS IS NO COORDINATES FIELD');
            if isempty(filename_data)
                filename_data = writepxdmf();
                filename_data.xdmf=1;
            end
            if ischar(filename_data)
                data = writepxdmf();
                data.xdmf = 1;
                data.filename = filename_data;
            else
                data = filename_data;
            end
            next_dim = numel(data.nodes)+1;
            
            if next_dim>1
                data.xdmf = 0;
            end
            mask = cellfun(@ischar,varargin);
            %options = varargin(mask);
            fields = [varargin{~mask}];
            
            %for i=1:numel(fields)
            %fields(i).set_current_domain(name);
            %end
            set_current_domain(obj,name);
            
            assert(obj.CURRENT_ID>0,'UNKNOWN DOMAIN FOR COORDINATE SYSTEM %s',obj.name);
            
            MAPPING = obj.get_current_elements();
            if ~isempty(fields)
            fields.set_current_domain(name);
            
            assert(all([obj.CURRENT_ID]>0),'UNKNOWN DOMAIN FOR SOME FIELDS');
            
            
            ELEMENTS = fields.get_current_elements();
            cell_field_mask = false(1,numel(fields));
            for i=1:numel(fields)
                cell_field_mask(i) = (numel(ELEMENTS(i).dof_entities)==1) && all(ELEMENTS(i).dof_entities==0);
            end
            
            if any(~cell_field_mask)
                if (next_dim==1)
                    data.nodes_fields = cell(1,nnz(~cell_field_mask));
                    data.nodes_fields_names = {fields(~cell_field_mask).name};
                else
                    assert(all(strcmp(data.nodes_fields_names,{fields(~cell_field_mask).name})),'Inconsistent Naming of node fields between dimensions');
                end
                data.nodes_fields{next_dim,nnz(~cell_field_mask)} = [];
                nxt = 0;
                for i=find(~cell_field_mask)
                    nxt = nxt+1;
                    if (fields(i).type==1) && strcmp(MAPPING.name,ELEMENTS(i).name)
                        tmp = fields(i).get_values(name);
                    else
                        if fields(i).type==1
                            tmp = obj.CURRENT_INTERP.interpolate(fields(i).CURRENT_INTERP,fields(i).get_values(name));
                        else
                            tmp = fields(i).FCT(obj.CURRENT_INTERP.interpolate(fields(i).PARENT_FIELD.CURRENT_INTERP,fields(i).PARENT_FIELD.get_values(name)));
                        end
                    end
                    if fields(i).Ncomp==1
                        data.nodes_fields{next_dim,nxt} = tmp';
                    else
                        for m = 1:fields(i).Nmodes
                            buf = tmp(:,(1:fields(i).Ncomp)+(m-1)*fields(i).Ncomp);
                            if size(buf,2)==2
                                buf = [buf zeros(size(buf,1),1)]; %#ok<AGROW>
                            end
                            data.nodes_fields{next_dim,nxt}(m,:) = reshape(buf',1,[]);
                        end
                    end
                end
            end
            
            if any(cell_field_mask)
                if (next_dim==1)
                    data.elements_fields = cell(1,nnz(cell_field_mask));
                    data.cell_fields_names = {fields(cell_field_mask).name};
                else
                    assert(all(strcmp(data.cell_fields_names,{fields(cell_field_mask).name})),'Inconsistent Naming of cell fields between dimensions');
                end
                data.elements_fields{next_dim,nnz(cell_field_mask)} = [];
                nxt = 0;
                for i=find(cell_field_mask)
                    nxt = nxt+1;
                    if fields(i).type == 1
                        tmp = fields(i).get_values(name);
                    else
                        tmp = fields(i).FCT(fields(i).PARENT_FIELD.get_values(name));
                    end
                    if fields(i).Ncomp==1
                        data.cell_fields{next_dim,nxt} = tmp';
                    else
                        for m = 1:fields(i).Nmodes
                            buf = tmp(:,(1:fields(i).Ncomp)+(m-1)*fields(i).Ncomp);
                            if size(buf,2)==2
                                buf = [buf zeros(size(buf,1))]; %#ok<AGROW>
                            end
                            data.cell_fields{next_dim,nxt}(m,:) = reshape(buf',1,[]);
                        end
                    end
                end
                
                
            end
            end
            %options
            X = obj.get_values(name);
            if size(X,2) > 1
                % the number of characters is the same as the number of dims
                if size(X,2) == size(obj.name,2) 
                    data.names{next_dim,1} = num2cell(obj.name);
                else
                    % if we have a name with the correct number of
                    % separators
                    if length(find(obj.name=='-')) == size(X,2)-1
                        data.names{next_dim,1} = regexp(obj.name, '-', 'split');
                    else
                        % in this case we just add a number to the name
                        data.names{next_dim,1} = strcat(obj.name,cellfun(@(arg) num2str(arg),num2cell(1:size(X,2)),'uniformoutput',false));
                    end
                end
            else
                 data.names{next_dim,1} = {obj.name};
            end
            
            if size(X,2) < 3
                X = [X zeros(size(X,1),3-size(X,2))];
            end
            data.nodes{next_dim,1} = X;
            
            
            data.from1 = 0;
            data.cells{next_dim,1} = renumber_compact(obj.CURRENT_INTERP.CONNECTIVITY)-1;
            data.cell_names{next_dim} = MAPPING.xdmf_name;
            if ischar(filename_data)
                writexdmf(data);
            end
        end
        
        function result = find_field(obj,field_name)
            %identifies a field by its name in a vector of fields
            result = find(strcmp(field_name,{obj.name}));
            assert(~isempty(result),'Could not find field');
        end
        
        function [field_id,dom_id,cols,comp_id,mode_id] = get_ids(obj,opts)
            %process the possible inputs of some function
            %identifies a field number, a domain number modes and
            %components from the input
            
            mask = cellfun(@isnumeric,opts);
            num = {};
            others = {};
            if any(mask)
                num = opts(mask);
            end
            if any(~mask)
                others = opts(~mask);
            end
            
            switch numel(others)
                case 0
                    error('You must specify a domain name')
                case 1
                    assert(numel(obj)==1,'Field name can only be ommited if there is one field')
                    field_id = 1;
                    dom_id = obj(field_id).INTERP.find_domain(others{1});
                case 2
                    field_id = obj.find_field(others{1});
                    dom_id = obj(field_id).INTERP.find_domain(others{2});
            end
            
            switch numel(num)
                case 0
                    comp_id = 1:obj(field_id).Ncomp;
                    mode_id = 1:obj(field_id).Nmodes;
                case 1
                    comp_id = num{1};
                    mode_id = 1:obj(field_id).Nmodes;
                case 2
                    comp_id = num{1};
                    mode_id = num{2};
            end
            assert(all(comp_id<=obj(field_id).Ncomp),'Error: request for a field component greater that the number of components');
            cols = reshape(bsxfun(@plus,comp_id',(mode_id-1)*obj(field_id).Ncomp),1,[]);
        end
        
        
    end
end