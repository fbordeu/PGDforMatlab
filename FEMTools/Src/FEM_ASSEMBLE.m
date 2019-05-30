%FEM_ASSEMBLE
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Adrien Leygue (Adrien.Leygue@ec-nantes.fr)
%
function varargout = FEM_ASSEMBLE(DOM_NAME,COORDS,U_FIELDS,T_FIELDS,K_FIELDS,INTEG_PARAM,varargin)
    WEAK = [];
    for i=1:numel(varargin)
        if ~ischar(varargin{i}), continue; end;
        switch upper(varargin{i})
            case 'LHS'
                WEAK = [WEAK varargin{i+1}]; %#ok<AGROW>
                continue;
            case 'RHS'
                WEAK = [WEAK varargin{i+1}]; %#ok<AGROW>
                continue;
            otherwise
                error(['Unrecognized option:' varargin{i}]);
        end
    end
    
    %basic check to speed up the assembly
    sym_galerkin = all(U_FIELDS==T_FIELDS);
    
    %gather all fields and check if different fields (handles) have
    %different names
    ALL_FIELDS = unique([COORDS; T_FIELDS(:);U_FIELDS(:);K_FIELDS(:)]);
    assert(numel(unique({ALL_FIELDS.name}))==numel(ALL_FIELDS),'FEM_ASSEMBLE: Different fields share the same name');
    
    %Prepare the integration over all the fields
    Nelem = ALL_FIELDS.set_current_domain(DOM_NAME);
    INTEG_WEIGHTS = ALL_FIELDS.prepare_integration(INTEG_PARAM);
    
    %compute the size of the local and global matrices
    %compute the assembly indexes in the local matrix
    [Nglob_col,Nloc_col] = U_FIELDS.compute_offsets;
    loc_cols = cell(1,numel(U_FIELDS));
    for i=1:numel(U_FIELDS)
        loc_cols{i} = U_FIELDS(i).get_local_index;
    end
    
    %sym_galerkin could be used here,... but it would not help much
    [Nglob_row,Nloc_row] = T_FIELDS.compute_offsets;
    loc_lines = cell(1,numel(T_FIELDS));
    for i=1:numel(T_FIELDS)
        loc_lines{i} = T_FIELDS(i).get_local_index;
    end
    
    %process the weak form to:
    % 1-extract LHS and RHS
    % 2-extract a unique element list to compute shape function gradients
    % 3-extract a catalog of elementary expressions that will be evaluated only once
    [LHS,RHS,ALL_ELEMENTS,catalog_phi,catalog_value] = process_weak_form(WEAK,U_FIELDS,T_FIELDS,ALL_FIELDS,loc_lines,loc_cols,Nloc_row,Nloc_col);
    
    %memory allocation for the assembly step
    has_rhs = ~isempty(RHS);
    has_lhs = ~isempty(LHS);
    if has_lhs
        sz_loc = Nloc_row*Nloc_col;
        LHS_i = zeros(Nelem*sz_loc,1);
        LHS_j = zeros(Nelem*sz_loc,1);
        LHS_values = zeros(Nelem*sz_loc,1);
    end
    if has_rhs
        RHS_i = zeros(Nelem*Nloc_row,1);
        RHS_values = zeros(Nelem*Nloc_row,1);
    end
    
    N_cat_phi = size(catalog_phi,1);
    N_cat_value = size(catalog_value,1);
    
    eval_catalog_phi = cell(1,N_cat_phi);
    eval_catalog_value = zeros(N_cat_value,numel(INTEG_WEIGHTS));
    
    rep_row = ones(1, Nloc_row);
    rep_col = ones(1, Nloc_col);
    
    %Loop over the elements
    for el=1:Nelem
        
        %local to global assembly map
        map_i = T_FIELDS.get_local_to_global_map(el);
        
        %should check the use of iJ
        [J,detJ] = COORDS.jacobian_integration(el);
        WJ = INTEG_WEIGHTS .* detJ';
        
        %iJ = cellfun(@inv,J,'uniformoutput',0);
        
        for i=1:N_cat_phi
            eval_catalog_phi{i} = ALL_ELEMENTS(catalog_phi(i,1)).gen_eval_phi(J,catalog_phi(i,2:end));
        end
        for i=1:N_cat_value
            eval_catalog_value(i,:) = ALL_FIELDS(catalog_value(i,1)).eval_field_integration(el,J,catalog_value(i,4:end),catalog_value(i,2));
        end
        
        if has_lhs
            A_loc = zeros(Nloc_row*Nloc_col,1);
            
            %local to global assembly map
            if sym_galerkin
                map_j = map_i;
            else
                map_j = U_FIELDS.get_local_to_global_map(el);
            end
            
            for tt = 1:numel(LHS)
                
                D = WJ*LHS(tt).factor; %scalar value
                %aux_fields
                if LHS(tt).has_known
                    D = D.* prod(eval_catalog_value(LHS(tt).known_entries,:),1);
                end
                
                tmp = eval_catalog_phi{LHS(tt).test_entry}*diag(D)*eval_catalog_phi{LHS(tt).unknown_entry}';
                A_loc(LHS(tt).loc_idx)=A_loc(LHS(tt).loc_idx) + tmp(:);
                
                
            end
            
            tmp = (1:sz_loc)+(el-1)*sz_loc;
            LHS_i(tmp) = map_i(rep_col,:)';
            LHS_j(tmp) = map_j(rep_row,:);
            LHS_values(tmp) = A_loc(:);
        end
        
        if has_rhs
            B_loc = zeros(Nloc_row,1);
            for tt = 1:numel(RHS)
                D = WJ*RHS(tt).factor; %scalar value
                %aux_fields
                if RHS(tt).has_known
                    D = D.* prod(eval_catalog_value(RHS(tt).known_entries,:),1);
                end
                B_loc(RHS(tt).loc_lines)=B_loc(RHS(tt).loc_lines) + eval_catalog_phi{RHS(tt).test_entry}*D';            
            end
            tmp = (1:Nloc_row)+(el-1)*Nloc_row;
            RHS_i(tmp) = map_i(:);
            RHS_values(tmp) = B_loc(:);
        end
        
    end
    
    varargout = cell(1,2);
    varargout{1} = sparse(Nglob_row,Nglob_col);
    varargout{2} = sparse(Nglob_row,1);
    if has_lhs
        varargout{1} = sparse(LHS_i,LHS_j,LHS_values,Nglob_row,Nglob_col);
    end
    if has_rhs
        varargout{2} = accumarray(RHS_i,RHS_values,[Nglob_row 1]);
    end 
    if nargout==1
        varargout = varargout([has_lhs has_rhs]);
    end
    
end

