function [LHS,RHS,ALL_ELEMENTS,catalog_phi,catalog_value] = process_weak_form(WEAK,U_FIELDS,T_FIELDS,ALL_FIELDS,loc_lines,loc_cols,Nloc_row,Nloc_col)
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Adrien Leygue (Adrien.Leygue@ec-nantes.fr)
%

    %remove terms with 0 as factor
    mask = ([WEAK.factor]~=0);
    WEAK = WEAK(mask);
    
    %extract all the field names for:
    U_names = {U_FIELDS.name}; %unknown fields
    T_names = {T_FIELDS.name}; %test fields
    ALL_names = {ALL_FIELDS.name}; %all fields
    
    %creation of a unique list of all elements
    tmp =ALL_FIELDS.get_current_elements(); %extract all elements
    ELEMENT_names = {tmp.name}; %get their names
    [~,gather_elts,scatter_elts] = unique(ELEMENT_names); %find the position of the unique names
    ALL_ELEMENTS = tmp(gather_elts); %gather the elements
    %scatter_elts will be used later to connect an alamant in ALL_FIELDS
    %with a particular element.
    
    %initialization of the structures for LHS & RHS
    TERMS = struct('factor',{},...
        'test_entry',{},'test_pos',{},'test_comp',{},...
        'unknown_entry',{},'unknown_pos',{},'unknown_comp',{},...
        'known_entries',{},'has_known',{},...
        'loc_lines',{},'loc_cols',{},'loc_idx',{});
    
    %initialization of the catalogs
    %catalogs regroup all the expressions that have to be evaluated
    %the catalog are first constructed with possible redundancy, unique
    %elements (rows) will be extracted later.
    catalog_phi = zeros(0,1+numel(WEAK(1).test_diff_symbol)); % element_position_in ALL_ELEMENTS diff_symbol
    catalog_value = zeros(0,1+1+1+numel(WEAK(1).test_diff_symbol)); % field_position_in_ALL_FIELDS component_number  mode_number diff_symbol
    idx_phi = 0; %index of the current catalog_phi entry
    idx_val = 0; %index of the current catalog_value entry
    for tt = 1:numel(WEAK) %loop over the terms
        
        TERMS(tt).factor = WEAK(tt).factor;
        
        %check if the test field exists
        [found,pos] = ismember(WEAK(tt).test_id,T_names);
        assert(found,['could not find ' WEAK(tt).test_id ') among the provided test fields']);
        TERMS(tt).test_pos = pos;
        TERMS(tt).test_comp = WEAK(tt).test_comp;
        TERMS(tt).loc_lines = loc_lines{pos}(WEAK(tt).test_comp,:);
        %add an entry to catalog_phi
        idx_phi = idx_phi+1;
        [~,pos] = ismember(WEAK(tt).test_id,ALL_names);
        entry_phi = [scatter_elts(pos) WEAK(tt).test_diff_symbol'];
        catalog_phi(idx_phi,:) = entry_phi;
        TERMS(tt).test_entry = idx_phi;
        
        if(~isempty(WEAK(tt).unknown_id))
            %check if the unknown field exists
            [found,pos] = ismember(WEAK(tt).unknown_id,U_names);
            assert(found,['could not find ' WEAK(tt).unknown_id ') among the provided unknown fields']);
            TERMS(tt).unknown_pos = pos;
            TERMS(tt).unknown_comp = WEAK(tt).unknown_comp;
            TERMS(tt).loc_cols = loc_cols{pos}(WEAK(tt).unknown_comp,:);
            %add an entry to catalog_phi
            idx_phi = idx_phi+1;
            [~,pos] = ismember(WEAK(tt).unknown_id,ALL_names);
            entry_phi = [scatter_elts(pos) WEAK(tt).unknown_diff_symbol'];
            catalog_phi(idx_phi,:) = entry_phi;
            TERMS(tt).unknown_entry = idx_phi;
            
            A = zeros(Nloc_row,Nloc_col);
            A(TERMS(tt).loc_lines,TERMS(tt).loc_cols) = 1;
            tmp = find(A);
            TERMS(tt).loc_idx = tmp(:);
            
        end
        
        %add the known functions
        if(~isempty(WEAK(tt).known_id))
            %check for existence
            [found,pos] = ismember(WEAK(tt).known_id,ALL_names);
            assert(all(found),'could not find some known field among the provided known fields');
            TERMS(tt).has_known = true;
            
            %add entries in catalog_value field_pod component mode(here 1) diff_symbol
            entries_value = [pos' (WEAK(tt).known_comp)' ones(size((WEAK(tt).known_comp)')) (WEAK(tt).known_diff_symbol)'];
            TERMS(tt).known_entries = idx_val+(1:numel(pos));
            catalog_value(idx_val+(1:numel(pos)),:) = entries_value;
            idx_val = idx_val+numel(pos);
        end
    end
    %extract unique elements of the catalogs and modifies the indexes in LHS and RHS
    [catalog_phi,~,scatter_phi] = unique(catalog_phi,'rows');
    [catalog_value,~,scatter_value] = unique(catalog_value,'rows');
    RHS_mask = false(1,numel(TERMS));
    for tt = 1:numel(TERMS)
        RHS_mask(tt) = isempty(TERMS(tt).unknown_entry);
        TERMS(tt).test_entry = scatter_phi(TERMS(tt).test_entry);
        TERMS(tt).unknown_entry = scatter_phi(TERMS(tt).unknown_entry);
        TERMS(tt).known_entries = scatter_value(TERMS(tt).known_entries);
    end
    RHS = TERMS(RHS_mask);
    LHS = TERMS(~RHS_mask);
end

% function [LHS,RHS,ALL_ELEMENTS,catalog_phi,catalog_value] = process_weak_form(WEAK,U_FIELDS,T_FIELDS,ALL_FIELDS,loc_lines,loc_cols,Nloc_row,Nloc_col)
%
%     %remove terms with 0 as factor
%     mask = ([WEAK.factor]~=0);
%     WEAK = WEAK(mask);
%
%     %extract all the field names for:
%     U_names = {U_FIELDS.name}; %unknown fields
%     T_names = {T_FIELDS.name}; %test fields
%     ALL_names = {ALL_FIELDS.name}; %all fields
%
%     %creation of a unique list of all elements
%     tmp =ALL_FIELDS.get_current_elements(); %extract all elements
%     ELEMENT_names = {tmp.name}; %get their names
%     [~,gather_elts,scatter_elts] = unique(ELEMENT_names); %find the position of the unique names
%     ALL_ELEMENTS = tmp(gather_elts); %gather the elements
%     %scatter_elts will be used later to connect an alamant in ALL_FIELDS
%     %with a particular element.
%
%     %initialization of the structures for LHS & RHS
%     LHS = struct('factor',{},'test_entry',{},'test_pos',{},'unknown_entry',{},'unknown_pos',{},'known_entries',{},'has_known',{},'loc_lines',{},'loc_cols',{},'loc_idx',{});
%     RHS = struct('factor',{},'test_entry',{},'test_pos',{},'known_entries',{},'has_known',{},'loc_lines',{});
%
%     %initialization of the catalogs
%     %catalogs regroup all the expressions that have to be evaluated
%     %the catalog are first constructed with possible redundancy, unique
%     %elements (rows) will be extracted later.
%     catalog_phi = zeros(0,1+numel(WEAK(1).test_diff_symbol)); % element_position_in ALL_ELEMENTS diff_symbol
%     catalog_value = zeros(0,1+1+numel(WEAK(1).test_diff_symbol)); % field_position_in_ALL_FIELDS component_number diff_symbol
%
%     idx_phi = 0; %index of the current catalog_phi entry
%     idx_val = 0; %index of the current catalog_value entry
%     idx_rhs = 0; %index of the current rhs element
%     idx_lhs = 0; %index of the current lhs element
%
%     for t = 1:numel(WEAK) %loop over the terms
%
%         if isempty(WEAK(t).unknown_id) %no unknown field --> rhs
%             %add an element in RHS
%             idx_rhs = idx_rhs+1;
%             RHS(idx_rhs).factor = WEAK(t).factor;
%
%             %check if the test field exists
%             [found,pos] = ismember(WEAK(t).test_id,T_names);
%             assert(found,['could not find ' WEAK(t).test_id ') among the provided test fields']);
%             RHS(idx_rhs).test_pos = pos;
%             RHS(idx_rhs).loc_lines = loc_lines{pos}(WEAK(t).test_comp,:);
%
%
%             %add an entry in catalog_phi
%             idx_phi = idx_phi+1;
%             [~,pos] = ismember(WEAK(t).test_id,ALL_names);
%             entry_phi = [scatter_elts(pos) WEAK(t).test_diff_symbol'];
%             catalog_phi(idx_phi,:) = entry_phi;
%             RHS(idx_rhs).test_entry = idx_phi;
%
%
%             %add the known functions
%             if(~isempty(WEAK(t).known_id))
%                 %check for existence
%                 [found,pos] = ismember(WEAK(t).known_id,ALL_names);
%                 assert(all(found),'could not find some known field among the provided unknown fields');
%                 RHS(idx_rhs).has_known = true;
%
%                 %add entries in catlog_value
%                 entries_value = [pos (WEAK(t).known_comp)' (WEAK(t).known_diff_symbol)'];
%                 RHS(idx_rhs).known_entries = idx_val+(1:numel(pos));
%                 catalog_value(idx_val+(1:numel(pos)),:) = entries_value;
%                 idx_val = idx_val+numel(pos);
%             end
%
%         else %we have a lhs term
%             %add an element in LHS
%             idx_lhs = idx_lhs+1;
%             LHS(idx_lhs).factor = WEAK(t).factor;
%
%             %check if the test field exists
%             [found,pos] = ismember(WEAK(t).test_id,T_names);
%             assert(found,['could not find ' WEAK(t).test_id ') among the provided test fields']);
%             LHS(idx_lhs).test_pos = pos;
%             LHS(idx_lhs).loc_lines = loc_lines{pos}(WEAK(t).test_comp,:);
%
%             %add an entry to catalog_phi
%             idx_phi = idx_phi+1;
%             [~,pos] = ismember(WEAK(t).test_id,ALL_names);
%             entry_phi = [scatter_elts(pos) WEAK(t).test_diff_symbol'];
%             catalog_phi(idx_phi,:) = entry_phi;
%             LHS(idx_lhs).test_entry = idx_phi;
%
%             %check if the unknown field exists
%             [found,pos] = ismember(WEAK(t).unknown_id,U_names);
%             assert(found,['could not find ' WEAK(t).unknown_id ') among the provided unknown fields']);
%             LHS(idx_lhs).unknown_pos = pos;
%             LHS(idx_lhs).loc_cols = loc_cols{pos}(WEAK(t).unknown_comp,:);
%
%             %add an entry to catalog_phi
%             idx_phi = idx_phi+1;
%             [~,pos] = ismember(WEAK(t).unknown_id,ALL_names);
%             entry_phi = [scatter_elts(pos) WEAK(t).unknown_diff_symbol'];
%             catalog_phi(idx_phi,:) = entry_phi;
%             LHS(idx_lhs).unknown_entry = idx_phi;
%
%             A = zeros(Nloc_row,Nloc_col);
%             A(LHS(idx_lhs).loc_lines,LHS(idx_lhs).loc_cols) = 1;
%             LHS(idx_lhs).loc_idx = find(A);
%
%             %add the known functions
%             if(~isempty(WEAK(t).known_id))
%                 %check for existence
%                 [found,pos] = ismember(WEAK(t).known_id,ALL_names);
%                 assert(all(found),'could not find some known field among the provided unknown fields');
%                 LHS(idx_lhs).has_known = true;
%
%                 %add entries in catalog_value
%                 entries_value = [pos (WEAK(t).known_comp)' (WEAK(t).known_diff_symbol)'];
%                 LHS(idx_lhs).known_entries = idx_val+(1:numel(pos));
%                 catalog_value(idx_val+(1:numel(pos)),:) = entries_value;
%                 idx_val = idx_val+numel(pos);
%             end
%         end
%     end
%
%     %extract unique elements of the catalogs and modifies the indexes in LHS and RHS
%     [catalog_phi,~,scatter_phi] = unique(catalog_phi,'rows');
%     [catalog_value,~,scatter_value] = unique(catalog_value,'rows');
%     for t = 1:numel(LHS)
%         LHS(t).test_entry = scatter_phi(LHS(t).test_entry);
%         LHS(t).unknown_entry = scatter_phi(LHS(t).unknown_entry);
%         LHS(t).known_entries = scatter_value(LHS(t).known_entries);
%     end
%     for t = 1:numel(RHS)
%         RHS(t).test_entry = scatter_phi(RHS(t).test_entry);
%         RHS(t).known_entries = scatter_value(RHS(t).known_entries);
%     end
% end
