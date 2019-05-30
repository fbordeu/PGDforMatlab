%FEM_ASSEMBLE_PGD
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Adrien Leygue (Adrien.Leygue@ec-nantes.fr)
%
function varargout = FEM_ASSEMBLE_PGD_2(DOM_NAME,COORDS,U_FIELDS,T_FIELDS,K_FIELDS,INTEG_PARAM,varargin)
    
    old_school = false;
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
            case 'OLD_SCHOOL'
                old_school = true;
                continue;
            otherwise
                error(['Unrecognized option:' varargin{i}]);
        end
    end
    %basic check to speed up the assembly
    %sym_galerkin = all(U_FIELDS==T_FIELDS);
    
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
    
    %TODO : ajouter les facteurs et faire l'expansion des aux-fields
    [LHS,RHS,ALL_ELEMENTS,catalog_phi,catalog_value] = process_weak_form(WEAK,U_FIELDS,T_FIELDS,ALL_FIELDS,loc_lines,loc_cols,Nloc_row,Nloc_col);
    has_rhs = ~isempty(RHS);
    has_lhs = ~isempty(LHS); 
    
    %virtual assembly and memory allocation
    if has_lhs
        size_term = zeros(1,numel(LHS));
        for tt = 1:numel(LHS)
            myterm = LHS(tt);
            my_i = myterm.test_entry; %test field number
            my_j = myterm.unknown_entry; %unnown field number
            size_term(tt) = ALL_ELEMENTS(catalog_phi(my_i,1)).Nphi * ALL_ELEMENTS(catalog_phi(my_j,1)).Nphi;
        end
        
        %memory allocation for virtual assembly
        LHS_i = zeros(Nelem*sum(size_term),1);
        LHS_j = zeros(Nelem*sum(size_term),1);
       
        idx_lhs = 0;
        for el = 1:Nelem
            for tt = 1:numel(LHS)
                myterm = LHS(tt);
                    
                map_i = U_FIELDS(myterm.test_pos).get_detailled_map(el,myterm.test_comp);
                map_j = T_FIELDS(myterm.unknown_pos).get_detailled_map(el,myterm.unknown_comp);
                
                tmp_i = zeros(numel(map_i),numel(map_j));
                tmp_j = zeros(numel(map_i),numel(map_j));
                
                for k = 1:numel(map_i)
                    tmp_j(k,:) = map_j;
                end
                for k = 1:numel(map_j)
                    tmp_i(:,k) = map_i;
                end
                tmp = idx_lhs+(1:(numel(map_i)*numel(map_j)));
                
                LHS_i(tmp) = map_i(ones(1,numel(map_j)),:)';
                LHS_j(tmp) = map_j(ones(1,numel(map_i)),:);
                idx_lhs = idx_lhs + numel(tmp);
            end
        end
        A = sparse(LHS_i,LHS_j,ones(size(LHS_i)),Nglob_row,Nglob_col);
        [pattern_i,pattern_j] = find(A);
        A_pattern = sparse(pattern_i,pattern_j,1:numel(pattern_i),Nglob_row,Nglob_col);
    end
    
    [LHS,RHS,catalog_value] = expand_separated_terms([LHS(:);RHS(:)],catalog_value,ALL_FIELDS);
    
    if has_lhs
        size_term = zeros(1,numel(LHS));
        for tt = 1:numel(LHS)
            myterm = LHS(tt);
            my_i = myterm.test_entry; %test field number
            my_j = myterm.unknown_entry; %unnown field number
            size_term(tt) = ALL_ELEMENTS(catalog_phi(my_i,1)).Nphi * ALL_ELEMENTS(catalog_phi(my_j,1)).Nphi;
        end 
        LHS_i = zeros(Nelem*sum(size_term),1);
        LHS_j = zeros(Nelem*sum(size_term),1);
        LHS_values = zeros(Nelem*sum(size_term),1);
    end
    
    if has_rhs
        size_term = zeros(1,numel(RHS));
        for tt = 1:numel(RHS)
            myterm = RHS(tt);
            my_i = myterm.test_pos;
            size_term(tt) = numel(myterm.loc_lines);%ALL_ELEMENTS(my_i).Nphi;
        end
        RHS_i = zeros(Nelem*sum(size_term),1);
        RHS_j = zeros(Nelem*sum(size_term),1);
        RHS_values = zeros(Nelem*sum(size_term),1);
    end

    N_cat_phi = size(catalog_phi,1);
    N_cat_value = size(catalog_value,1);
    
    eval_catalog_phi = cell(1,N_cat_phi);
    eval_catalog_value = zeros(N_cat_value,numel(INTEG_WEIGHTS));
    
    %Loop over the elements
    idx_lhs = 0;
    idx_rhs = 0;
     for el=1:Nelem
         
         [J,detJ] = COORDS.jacobian_integration(el);
         WJ = INTEG_WEIGHTS .* detJ';
         
         for i=1:N_cat_phi
             eval_catalog_phi{i} = ALL_ELEMENTS(catalog_phi(i,1)).gen_eval_phi(J,catalog_phi(i,2:end));
         end
         for i=1:N_cat_value
             eval_catalog_value(i,:) = ALL_FIELDS(catalog_value(i,1)).eval_field_integration(el,J,catalog_value(i,4:end),catalog_value(i,2),catalog_value(i,3));
         end
         
         if has_lhs
             for tt = 1:numel(LHS)               
               
                D = WJ*LHS(tt).factor; %scalar value
                %aux_fields
                if LHS(tt).has_known
                    D = D.* prod(eval_catalog_value(LHS(tt).known_entries,:),1);
                end
                
                contrib = eval_catalog_phi{LHS(tt).test_entry}*diag(D)*eval_catalog_phi{LHS(tt).unknown_entry}';
                
                map_i = T_FIELDS(LHS(tt).test_pos).get_detailled_map(el,LHS(tt).test_comp);
                map_j = U_FIELDS(LHS(tt).unknown_pos).get_detailled_map(el,LHS(tt).unknown_comp);
                
                tmp = A_pattern(map_i,map_j);
                indices = idx_lhs+(1:numel(map_i)*numel(map_j));
                
                LHS_i(indices) = tmp(:);
                LHS_j(indices) = tt;
                LHS_values(indices) = contrib(:);
                idx_lhs = idx_lhs + numel(indices);
            end
         end
         if has_rhs
             for tt = 1:numel(RHS)
                 D = WJ*RHS(tt).factor; %scalar value
                 %aux_fields
                 if RHS(tt).has_known
                     D = D.* prod(eval_catalog_value(RHS(tt).known_entries,:),1);
                 end
                 contrib =  eval_catalog_phi{RHS(tt).test_entry}*D';
                 
                 map_i = T_FIELDS(RHS(tt).test_pos).get_detailled_map(el,RHS(tt).test_comp);
                 
                 indices = idx_rhs+(1:numel(map_i));
                 RHS_i(indices) = map_i;
                 RHS_j(indices) = tt;
                 RHS_values(indices) = contrib(:);
                 idx_rhs = idx_rhs + numel(indices);
             end
         end
         
     end
    
    
    if old_school
        varargout  = cell(1,double(has_lhs)+double(has_rhs));
        offset = 0;
        if has_lhs
            [tmp_i,tmp_j] = find(A_pattern);
            AA = cell(1,numel(LHS));
            ALL_ops = sparse(LHS_i,LHS_j,LHS_values,nnz(A_pattern),numel(LHS));
            for i=1:numel(LHS)
                AA{i} = sparse(tmp_i,tmp_j,ALL_ops(:,i),size(A_pattern,1),size(A_pattern,2));
            end
            varargout{1} = AA;
            offset = 1;
        else
            if nargout==2
                offset = 1;
                varargout{1} = cell(1,0);
            end
        end
        
        if has_rhs
            varargout{offset+1}  = accumarray([RHS_i RHS_j],RHS_values,[Nglob_row,numel(RHS)]);
        elseif nargout>offset
            varargout{offset+1} = zeros(size(A_pattern,1),1);
        end
            
        
    else
        varargout  = cell(1,3*double(has_lhs)+double(has_rhs));
        offset = 0;
        if has_lhs
            [tmp_i,tmp_j] = find(A_pattern);
            varargout(1:4) = {sparse(LHS_i,LHS_j,LHS_values,nnz(A_pattern),numel(LHS)),tmp_i,tmp_j,size(A_pattern)};
            offset = 4;
        end
        
        if has_rhs
            varargout{offset+1}  = accumarray([RHS_i RHS_j],RHS_values,[Nglob_row,numel(RHS)]);
        end
    end    
end

