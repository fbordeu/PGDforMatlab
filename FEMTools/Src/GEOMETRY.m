%GEOMETRY CLASS
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Adrien Leygue (Adrien.Leygue@ec-nantes.fr)
%
%SYM CLASSES PACKAGE
%
%
%Last modification of this file: 23 sept. 2013
classdef GEOMETRY
    properties
        ELEMENTS
        ENTITIES
        SUB_ENTITIES
        PHYSICAL_TAG
        PHYSICAL_NAMES
        NUMBER_TAG
        DOF_NUMBERING
        %KEYS
        NELEM
        Default
    end
    methods
        function obj = GEOMETRY(TOPO_p,NUMBER_TAG_p,PHYSICAL_TAG_p,PHYSICAL_NAMES_p,ELEMENTS_NAMES_p)
            %One fills the structure with exiting data
            obj.ELEMENTS = ELEMENT(ELEMENTS_NAMES_p,true);
            Ntypes = numel(obj.ELEMENTS);
            
            [~,pos] = ismember(ELEMENTS_NAMES_p,{obj.ELEMENTS.name});
            
            obj.ENTITIES = cell(1,Ntypes);
            obj.NUMBER_TAG = cell(1,Ntypes);
            obj.PHYSICAL_TAG = cell(1,Ntypes);
            
            obj.ENTITIES(pos) = TOPO_p;
            obj.NUMBER_TAG(pos) = NUMBER_TAG_p;
            obj.PHYSICAL_TAG(pos) = PHYSICAL_TAG_p;
            
            obj.PHYSICAL_NAMES = PHYSICAL_NAMES_p;
            
            %obj.KEYS = cell(Ntypes,1);
            KEYS = cell(Ntypes,1);
            
            obj.DOF_NUMBERING = cell(Ntypes,1);
            obj.SUB_ENTITIES = cell(Ntypes,1);
            
            
            %Keys allow the unique identification of an entity (ordered set of nodes)
            %Once all entities are created, one preferably uses the DOF_NUMBERING quantity to uniquely identify an entity
            for i=1:Ntypes
                %obj.KEYS{i} = sort(obj.ENTITIES{i},2);
                KEYS{i} = make_key(obj.ENTITIES{i});
            end
            obj.NELEM = zeros(1,Ntypes);
            %Main loop to generate all the unique sub-entities
            %and start filling the DOF_NUMBERING field
            %elements are processed by decreasing dimensionality
            for elem = Ntypes:-1:1
                if size(obj.ENTITIES{elem},1)==0
                    continue;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%
                %First remove duplicates
                %%%%%%%%%%%%%%%%%%%%%%%%
                %identifies the unique keys while preserving the order, therefore this should not delete the original entities
                %[~,gather_vect] = unique(obj.KEYS{elem},'rows','stable');
                [~,gather_vect] = unique(KEYS{elem},'rows','stable');
                
                %one has to check that the discarded elements do not belong to the original elements set
                %as the NUMBER_TAG is added last, we check that the NUMBER_TAG of the discarded elements is 0
                %mask = true(size(obj.KEYS{elem},1),1);
                mask = true(size(KEYS{elem},1),1);
                
                mask(gather_vect) = false;
                tmp = obj.NUMBER_TAG{elem}(mask,:);
                assert(all(~tmp(:)),'Error, trying to remove previously existing elements');
                
                %remove duplicates
                obj.NUMBER_TAG{elem} = obj.NUMBER_TAG{elem}(gather_vect,:);
                obj.PHYSICAL_TAG{elem} = obj.PHYSICAL_TAG{elem}(gather_vect,:);
                obj.ENTITIES{elem} = obj.ENTITIES{elem}(gather_vect,:);
                %obj.KEYS{elem} = obj.KEYS{elem}(gather_vect,:);
                KEYS{elem} = KEYS{elem}(gather_vect,:);
                
                Nelem = size(obj.ENTITIES{elem},1);
                %%%%%%%%%%%%%%%%%%
                %Add elements generated by the current elements
                %%%%%%%%%%%%%%%%%%
                
                %create DOF_NUMBER for the entities
                obj.DOF_NUMBERING{elem} = zeros(Nelem,1);
                %check if this element has a DOF associated to itself and finds the position
                [~, dof_number] = ismember(0,obj.ELEMENTS(elem).dof_entities);
                %if ~isempty(dof_number)
                if dof_number
                    obj.DOF_NUMBERING{elem} = obj.ENTITIES{elem}(:,dof_number); % give the entity number to the DOF
                end
                %other dofs not associated with the element itself will be numbered later
                %Generate sub entities
                for sub_el=1:obj.ELEMENTS(elem).Nentities
                    %Generate sub entities
                    [sub_name,tmp] = entity_table(obj.ELEMENTS(elem),obj.ENTITIES{elem},sub_el);
                    %sub_type = obj.ELEMENTS.identify(sub_name);
                    sub_type = strcmpi(sub_name,{obj.ELEMENTS.name});
                    %add to existing entities
                    obj.ENTITIES{sub_type} = [obj.ENTITIES{sub_type}; tmp];
                    %obj.KEYS{sub_type} = [obj.KEYS{sub_type}; obj.make_key(tmp)];
                    KEYS{sub_type} = [KEYS{sub_type}; make_key(tmp)];
                    
                    %DOF_NUMBERING associated with the sub_entities will be generated when those
                    %sub_elements are processed
                    %TAGs will be generated later only memory allocation
                    obj.NUMBER_TAG{sub_type} = [obj.NUMBER_TAG{sub_type}; zeros(Nelem,1)];
                    obj.PHYSICAL_TAG{sub_type} = [obj.PHYSICAL_TAG{sub_type}; zeros(Nelem,1)];
                end
            end
            
            %At this stage: all entities exist and have a correct self-referent DOF_NUMBERING
            %have a key + PHYSICAL_TAG (PHYSICAL_TAG is zero for new entities)
            %the NUMBER_TAG need to be added
            %for each element we still have to identify the DOF_NUMBERING
            %of its sub_entities
            
            %generate new NUMBERTAG and DOF_NUMBER--> continued numbering
            mask = ~cellfun(@isempty,obj.NUMBER_TAG);
            new_number_tag = max(cellfun(@max,obj.NUMBER_TAG(mask)))+1;
            mask = ~cellfun(@isempty,obj.DOF_NUMBERING);
            new_DOF_NUMBER = max(cellfun(@max,obj.DOF_NUMBERING(mask)))+1;
            
            %Process the elements in order of increasing dimensionality
            %TO CORRECT (origilal= 1:Ntypes)
            for elem=1:Ntypes %obj.ELEMENTS.sort_dim('ascend')
                if size(obj.ENTITIES{elem},1)==0
                    continue;
                end
                obj.Default = elem;
                %Identifies the unnumbered entities
                mask = obj.NUMBER_TAG{elem}==0;
                mask_dof = obj.DOF_NUMBERING{elem}==0;
                N_new = nnz(mask);
                N_new_dof = nnz(mask_dof);
                %Numbering of the unnumbered entities
                obj.NUMBER_TAG{elem}(mask) = new_number_tag:(new_number_tag+N_new-1);
                obj.DOF_NUMBERING{elem}(mask_dof) = new_DOF_NUMBER:(new_DOF_NUMBER+N_new_dof-1);
                %prepare for next numbering
                new_number_tag = new_number_tag+N_new;
                new_DOF_NUMBER = new_DOF_NUMBER+N_new_dof;
                %builds the table identifying sub-entities by their number
                Nelem = size(obj.ENTITIES{elem},1);
                obj.SUB_ENTITIES{elem} = zeros(Nelem,obj.ELEMENTS(elem).Nentities);
                %Loop over the sub-entities
                for sub_el=1:obj.ELEMENTS(elem).Nentities
                    
                    [sub_name,tmp] = entity_table(obj.ELEMENTS(elem),obj.ENTITIES{elem},sub_el);
                    %identify sub-entity type
                    %sub_type = obj.ELEMENTS.identify(sub_name);
                    sub_type = strcmpi(sub_name,{obj.ELEMENTS.name});
                    %matches the sub-entities to existing entities
                    %[~,pos] = ismember(obj.make_key(tmp),obj.KEYS{sub_type},'rows');
                    [~,pos] = ismember(make_key(tmp),KEYS{sub_type},'rows');
                    
                    %fill the corresponding column
                    obj.SUB_ENTITIES{elem}(:,sub_el) = obj.DOF_NUMBERING{sub_type}(pos);
                end
            end
            obj.NELEM = cellfun(@numel,obj.NUMBER_TAG);
            
        end
        
        function DOF_ENTITIES = fill_DOF_ENTITIES(GEOMETRY,mapping_id,local_entities,mask)
            %fill a connectivity table for an element with dofs described by
            %local_entities living on some mapping mesh + mask
            DOF_ENTITIES = zeros(nnz(mask),numel(local_entities));
            if numel(unique(local_entities)) == numel(local_entities)
                if any(local_entities==0)
                    DOF_ENTITIES(:,local_entities==0)= GEOMETRY.DOF_NUMBERING{mapping_id}(mask);
                end
                DOF_ENTITIES(:,local_entities~=0) = GEOMETRY.SUB_ENTITIES{mapping_id}(mask,local_entities(local_entities~=0));
            else
                for i=1:numel(local_entities)
                    if local_entities(i)==0
                        DOF_ENTITIES(:,i)= GEOMETRY.DOF_NUMBERING{mapping_id}(mask)+ nnz(local_entities(1:(i-1))==local_entities(i))*sum(GEOMETRY.NELEM);
                    else
                        DOF_ENTITIES(:,i) = GEOMETRY.SUB_ENTITIES{mapping_id}(mask,local_entities(i))+ nnz(local_entities(1:(i-1))==local_entities(i))*sum(GEOMETRY.NELEM);
                    end
                end
            end
        end
        
        
        
    end
end
function result = make_key(arg)
    result = sort(arg,2);
end