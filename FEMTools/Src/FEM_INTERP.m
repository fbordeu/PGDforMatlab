%FEM_INTERP CLASS
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Adrien Leygue (Adrien.Leygue@ec-nantes.fr)
%
classdef FEM_INTERP < handle
    properties
        ELEM = ELEMENT;
        DOMAIN = '';
        %PARENT = [];
        N_ELEMENTS = 0;
        N_VALUES = 0;
        MAX_VALUES = 0;
        CONNECTIVITY = [];
        NUMBER_TAG = [];
        PHYSICAL_TAG = [];
        
        NPCIntegrationPoints = [];
        PCIntegrationParam = [];
        NPCIntegrationPointsPerElement = [];
        PCIntegrationWeights = {};
        
        PCMatrix = {};
        PCTableValues = {};
        PCTableIndex = {};
        
        dofs = [];
        initdofs = false;
    end
    
    methods
        function obj = FEM_INTERP(GEOMETRY,DOMAIN,EL_NAME,options)
            %domain: un domaine ou bien ALL
            %Option SUB seulement si ALL
            
            if nargin==0 , return; end
            %Handle possible empty or missing argument
            if nargin<2
                DOMAIN='ALL';
            end
            if nargin<3
                EL_NAME = 'DEFAULT';
            end
%             if isempty(DOMAIN)
%                 DOMAIN='ALL';
%             end
%             if isempty(EL_NAME)
%                 EL_NAME = 'DEFAULT';
%             end
            
            %Validate inputs
            assert(isa(GEOMETRY,'GEOMETRY'),'First argument has to be of type GEOMETRY');
            assert(ischar(DOMAIN),'Second argument has to be a string (or empty)');
            assert(ischar(EL_NAME),'Second argument has to be a string (or empty)');           
            
            %Parse options
            discontinuous=0;
            if nargin==4
                if ~iscell(options)
                    options = {options};
                end
                assert(all(cellfun(@ischar,options)),'Fourth argument has to be a string or a cell array of strings')
                discontinuous = any(strcmpi('Discontinuous',options));
            else
                options = {};
            end
            
            %identify domain, mapping element, elements of the mesh
            if  any(strcmpi(DOMAIN,'ALL')) %If 'ALL' is used for the domain
                DOMAIN = 'ALL';
                mapping_id = GEOMETRY.Default;
                ELEMENT_MAPPING = GEOMETRY.ELEMENTS(mapping_id);
                mask = true(GEOMETRY.NELEM(mapping_id),1);
            else %otherwise if a specific domain is given
                my_physical_tag = find(strcmpi(DOMAIN, GEOMETRY.PHYSICAL_NAMES));
                assert(~isempty(my_physical_tag),'Unknown DOMAIN');
                %find in which element types the domain is (should be one only)
                idx = find(cellfun(@(in)  any(in==my_physical_tag),GEOMETRY.PHYSICAL_TAG));
                assert(~isempty(idx),'Empty Domain')
                assert(numel(idx)==1,'Domain contains different element types... not supported yet');
                %id the mapping
                mapping_id = idx;
                ELEMENT_MAPPING = GEOMETRY.ELEMENTS(mapping_id);
                %select elements belonging to the domain
                mask = GEOMETRY.PHYSICAL_TAG{mapping_id}==my_physical_tag;               
            end
            %id the interpolation
            if ~strcmpi(EL_NAME,'DEFAULT')
                ELEMENT_INTERP = ELEMENT(EL_NAME);
                assert(check_shape(ELEMENT_INTERP,ELEMENT_MAPPING),'Incompatible elements shapes')
            else
                ELEMENT_INTERP = ELEMENT_MAPPING;
            end
            %Fill the connectivity of the domain
            if(discontinuous) %in case of discontinuous interpolation, each element can be numbered independently
                obj.N_VALUES = ELEMENT_INTERP.Nphi *nnz(mask);
                DOM_CONNECTIVITY = reshape(1:obj.N_VALUES,[ELEMENT_INTERP.Nphi nnz(mask)])';
            else %continuous interpolation, dofs of different elements have to be matched to a single geometrical entity
                %Generation of the numbering,...we map each local dof with a
                %unique geometrical entity
                DOF_ENTITIES = GEOMETRY.fill_DOF_ENTITIES(mapping_id,ELEMENT_INTERP.dof_entities,mask);
                %then we identify duplicate geometrical entities to match dofs
                %accross elements
                [DOM_CONNECTIVITY, obj.N_VALUES] = renumber_compact(DOF_ENTITIES);
            end
            
            %Output
            obj.DOMAIN = DOMAIN;
            %obj.PARENT = [];
            obj.ELEM = ELEMENT_INTERP;
            obj.N_ELEMENTS = size(DOM_CONNECTIVITY,1);
            obj.NUMBER_TAG = GEOMETRY.NUMBER_TAG{mapping_id}(mask);
            obj.PHYSICAL_TAG = GEOMETRY.PHYSICAL_TAG{mapping_id}(mask);
            obj.CONNECTIVITY = DOM_CONNECTIVITY;
            obj.N_VALUES = numel(unique(DOM_CONNECTIVITY));
            obj.MAX_VALUES = obj.N_VALUES;

            %If no subdomains, exit
            pos_options = find(strcmpi('SUB',options));
            if isempty(pos_options), return; end
            
            %assert that there is something after each 'SUB' in the options
            assert(numel(options)>max(pos_options),'inconsistent options for the subdomain selection');
            %initialize array 
            obj(numel(pos_options)+1) = FEM_INTERP;
            
            %Loop over all 'SUB' options
            for sub = 1:numel(pos_options)
                SUB_DOMAIN = options{pos_options(sub)+1};
              
                %identify the sub-domain, sub-element (mapping and interp) & elements of the sub-mesh
                my_sub_tag = find(strcmpi(SUB_DOMAIN, GEOMETRY.PHYSICAL_NAMES));
                assert(~isempty(my_sub_tag),'Unknown SUBDOMAIN');
                %find in which element types the sub-domain is (should be one only)
                idx = find(cellfun(@(in)  any(in==my_sub_tag),GEOMETRY.PHYSICAL_TAG));
                assert(~isempty(idx),'Empty SubDomain');
                assert(numel(idx)==1,'Sub-Domains contain different element types');
                mapping_id_sub = idx;
                ELEMENT_MAPPING_SUB = GEOMETRY.ELEMENTS(mapping_id_sub);
                mask_sub = GEOMETRY.PHYSICAL_TAG{mapping_id_sub} ==my_sub_tag;
                %when the subdomain is of the same type as the domain,...
                %we just select the elements belonging to the subdomain
                if (mapping_id_sub==mapping_id)
                    %identify the elements of the sub_domain
                    assert(all(mask(mask_sub)),'The subdomain is does not entirely belong to the domain');   
                    obj(sub+1).DOMAIN = SUB_DOMAIN; %#ok<*AGROW>
                    %obj(sub+1).PARENT = obj(1);
                    obj(sub+1).ELEM = ELEMENT_INTERP;
                    obj(sub+1).N_ELEMENTS = nnz(mask_sub);
                    obj(sub+1).NUMBER_TAG = GEOMETRY.NUMBER_TAG{mapping_id}(mask_sub);
                    obj(sub+1).PHYSICAL_TAG = GEOMETRY.PHYSICAL_TAG{mapping_id}(mask_sub);
                    obj(sub+1).CONNECTIVITY = obj(1).CONNECTIVITY(mask_sub,:);
                    obj(sub+1).N_VALUES = numel(unique(obj(sub+1).CONNECTIVITY));
                    obj(sub+1).MAX_VALUES = obj(1).MAX_VALUES;
                    continue;
                end
                %IF NOT,... WE HAVE A BOUNDARY OF SOME SORT
                %We have to check if it is generated by the domain
                sub_entities_pos = ELEMENT_MAPPING.identify_entity(ELEMENT_MAPPING_SUB.name);
                generated_sub =  ismember(GEOMETRY.DOF_NUMBERING{mapping_id_sub},GEOMETRY.SUB_ENTITIES{mapping_id}(mask,sub_entities_pos));
                assert( all(generated_sub(mask_sub)),'The subdomain is not generated by the domain');
                
                %here we should check for continuity and maybe assign
                %another element type
                interp_name_sub = ELEMENT_INTERP.entity_name(sub_entities_pos(1));
                ELEMENT_INTERP_SUB = ELEMENT(interp_name_sub);
                N_EL = nnz(mask_sub);
                
                if(discontinuous)
                    %check if a subdomain element appears more than once,...
                    %this would mean a discontinuous interface
                    idx_mask_sub = find(mask_sub);
                    %counts the number of times an element of the subdomain
                    %is generated by the DOMAIN's elements
                    occurences = zeros(N_EL,1);
                    for i=1:N_EL
                        my_sub_entity = GEOMETRY.DOF_NUMBERING{mapping_id_sub}(idx_mask_sub(i));
                        occurences(i) = nnz(GEOMETRY.SUB_ENTITIES{mapping_id}(mask,:) == my_sub_entity);
                    end
                    %For the moment we allow only 1 occurence therwise the
                    %element would be an interface element
                    assert(max(occurences)==1,'Discontinuous interface description not implemented');
                    %Fill the connectivity table
                    SUB_CONNECTIVITY = zeros(N_EL,ELEMENT_INTERP_SUB.Nphi*max(occurences));
                    for i=1:N_EL
                        %identify the identif of the sub_element
                        my_sub_entity = GEOMETRY.DOF_NUMBERING{mapping_id_sub}(idx_mask_sub(i));
                        %find the mother element (row) and the local entity
                        %number in the parent element(col)
                        [row,col] = find(GEOMETRY.SUB_ENTITIES{mapping_id}==my_sub_entity);
                        %generates the connectivity of the sub-element
                        %SUB_CONNECTIVITY(i,:) = ELEMENT_INTERP_SUB.entity_table(DOM_CONNECTIVITY(row,:),col);
                        [~, SUB_CONNECTIVITY(i,:)] = ELEMENT_INTERP.entity_table(DOM_CONNECTIVITY(row,:),col);
                        
                    end
                else
                    %Generation of the numbering,...we map each local dof with a
                    %unique geometrical entity
                    SUB_DOF_ENTITIES = GEOMETRY.fill_DOF_ENTITIES(mapping_id_sub,ELEMENT_INTERP_SUB.dof_entities,mask_sub);                    
                    [~,gather_vect] = unique(DOM_CONNECTIVITY,'sorted');
                    val_dof = DOF_ENTITIES(gather_vect);
                    [~,SUB_CONNECTIVITY] = ismember(SUB_DOF_ENTITIES,val_dof);
                end
                
                obj(sub+1).DOMAIN = SUB_DOMAIN; %#ok<*AGROW>
                %obj(sub+1).PARENT = obj(1);
                obj(sub+1).ELEM = ELEMENT_INTERP_SUB;
                obj(sub+1).N_ELEMENTS = N_EL;
                obj(sub+1).NUMBER_TAG = GEOMETRY.NUMBER_TAG{mapping_id_sub}(mask_sub);
                obj(sub+1).PHYSICAL_TAG = GEOMETRY.PHYSICAL_TAG{mapping_id_sub}(mask_sub);
                obj(sub+1).CONNECTIVITY = SUB_CONNECTIVITY;
                obj(sub+1).N_VALUES = numel(unique(obj(sub+1).CONNECTIVITY));
                obj(sub+1).MAX_VALUES = obj(1).MAX_VALUES;
                
            end
                    
        end % FUNCTION
        
        function result = get_dofs(obj)
                    
            if obj.initdofs
                result = obj.dofs;
            else
                obj.initdofs = true;
                obj.dofs = unique(obj.CONNECTIVITY(:));
                result =  obj.dofs;
            end
        end
        function result = get_element_dofs(obj,el)
            result = obj.CONNECTIVITY(el,:);
        end
        function result = get_Nphi(el)
            if nargin==0
                el = 1;
            end
            result = numel(obj.connectivity(el,:));
        end
        function result = find_domain(obj,name)
            result = find(strcmpi(name,{obj.DOMAIN}));
            if isempty(result)
                result = 0;
            end
        end
        
        function result = interpolate(obj,COORDS_INTERP,X)
            result = zeros(obj.N_VALUES,size(X,2));
            M  = eval_shape_functions(COORDS_INTERP.ELEM,obj.ELEM.nodes);
            tmp1 = renumber_compact(COORDS_INTERP.CONNECTIVITY);
            tmp2 = renumber_compact(obj.CONNECTIVITY);
            for el = 1:obj.N_ELEMENTS
                %Xi = X(COORDS_INTERP.CONNECTIVITY(el,:),:);
                Xi = X(tmp1(el,:),:);
                
                XX = M*Xi;
                %result(obj.CONNECTIVITY(el,:),:) = XX;
                result(tmp2(el,:),:) = XX;
                
            end
        end
        
        function result = merge_domains(obj,name1,name2,new_name)
            
            pos1 = find_domain(obj,name1);
            pos2 = find_domain(obj,name2);
            assert(pos1>0,'Unknown first domain');
            assert(pos2>0,'Unknown second domain');
            
            assert(strcmp(obj(pos1).ELEM.name,obj(pos2).ELEM.name),'Error: the two domains have different element types');
            assert(~any(strcmp(new_name,{obj.DOMAIN})),'Error: already existing new domain name');
            
            result = obj;
            result(numel(result)+1) = FEM_INTERP;
            result(end).ELEM = result(pos1).ELEM;
            result(end).DOMAIN = new_name;
            %result(end).PARENT = result(pos1).PARENT;
            
            [~,ia] = setdiff(result(pos1).NUMBER_TAG,result(pos2).NUMBER_TAG);
            mask2 = true(result(pos2).N_ELEMENTS,1);
            mask1 = false(result(pos1).N_ELEMENTS,1);
            mask1(ia) = true;
            
            result(end).CONNECTIVITY = [result(pos1).CONNECTIVITY(mask1,:);result(pos2).CONNECTIVITY(mask2,:)];
            
            result(end).N_ELEMENTS = size(result(end).CONNECTIVITY,1);
            result(end).N_VALUES = numel(unique(result(end).CONNECTIVITY(:)));
            result(end).MAX_VALUES = result(pos1).MAX_VALUES;
            
            result(end).NUMBER_TAG = [result(pos1).NUMBER_TAG(mask1);result(pos2).NUMBER_TAG(mask2)];
            result(end).PHYSICAL_TAG = [result(pos1).PHYSICAL_TAG(mask1);result(pos2).PHYSICAL_TAG(mask2)];
        end

        function [result,new_X] = extract_domain(obj,X,master,varargin)
            N_sub = numel(varargin);
            if N_sub ; assert(all(cellfun(@ischar,varargin)),'All extra arguments must be strings describing the sub-domains'); end
            
            pos = find_domain(obj,master);
            assert(pos>0,'Could not find domain to extract');
            pos_sub = zeros(1,N_sub);
            for i=1:N_sub
                pos_sub(i) = find_domain(obj,varargin{i});
                assert(pos_sub(i)>0,['Could not find sub-domain ' varargin{i}]);
            end
            all_pos = [pos pos_sub];
            
            result(numel(all_pos)) = FEM_INTERP;
            
            old_dofs = unique(obj(all_pos(1)).CONNECTIVITY(:),'sorted');
            new_X = X(old_dofs,:);
            %new_dofs = renumber_compact(old_dofs);
            for i=1:numel(all_pos)
                result(i).ELEM = obj(all_pos(i)).ELEM;
                result(i).DOMAIN = obj(all_pos(i)).DOMAIN;
                result(i).DOMAIN
                %result(i).PARENT = result(1);
                result(i).N_ELEMENTS = obj(all_pos(i)).N_ELEMENTS;
                result(i).N_VALUES = obj(all_pos(i)).N_VALUES;
                result(i).MAX_VALUES = numel(old_dofs);
                [~,SUB_CONNECTIVITY] = ismember(obj(all_pos(i)).CONNECTIVITY,old_dofs);
                assert(all(all(SUB_CONNECTIVITY)),'A subdomain has nodes absent from the main domain');
                result(i).CONNECTIVITY = SUB_CONNECTIVITY;
                result(i).NUMBER_TAG = obj(all_pos(i)).NUMBER_TAG;
                result(i).PHYSICAL_TAG = obj(all_pos(i)).PHYSICAL_TAG;
            end   
        end
        
    end %METHODS
    
    methods (Static=true)
        function result = simple_constructor(TOPO,NUMBER_TAG,PHYSICAL_TAG,PHYSICAL_NAMES,elem_names)
            result(numel(PHYSICAL_NAMES)) = FEM_INTERP;
            tmp = [];
            for i=1:numel(elem_names)
                dofs = TOPO{i};
                dofs = unique(dofs(:));
                result(i).N_VALUES = numel(dofs);
                tmp = [tmp;dofs];
            end
            tmp = sort(unique(tmp));
            assert(all(tmp'==(1:numel(tmp))),'Error: non-compact FEM numbering');
            all_physical = all(cellfun(@(arg) all(ismember(arg,1:numel(PHYSICAL_NAMES))),PHYSICAL_TAG));
            assert(all_physical,'Error some elements do ot seem to belong to any physical domain');
            
            for i=1:numel(PHYSICAL_NAMES)
                pos = find(cellfun(@(arg) any(arg==i),PHYSICAL_TAG));
                assert(~isempty(pos),'at least one domain does not have any element');
                assert(numel(pos)==1,'at least one domain contains elements of different types');
                
                mask = PHYSICAL_TAG{pos}==i;
                
                result(i).ELEM = ELEMENT(elem_names{pos});
                result(i).DOMAIN = PHYSICAL_NAMES{i};
                %result(i).PARENT = [];
                
                result(i).N_ELEMENTS = nnz(mask);
                result(i).CONNECTIVITY = TOPO{pos}(mask,:);
                result(i).N_VALUES = numel(unique(result(i).CONNECTIVITY(:)));
                result(i).MAX_VALUES = numel(tmp);                
                result(i).NUMBER_TAG = NUMBER_TAG{pos}(mask);
                result(i).PHYSICAL_TAG = PHYSICAL_TAG{pos}(mask);
            end
        end
    end
    
end %CLASSDEF
