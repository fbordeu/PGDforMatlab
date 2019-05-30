function [LHS,RHS,catalog_value] = expand_separated_terms(WEAK,old_catalog_value,ALL_FIELDS)
%     result = struct('factor',{},...
%         'test_entry',{},'test_pos',{},'test_comp',{},...
%         'unknown_entry',{},'unknown_pos',{},'unknown_comp',{},...
%         'known_entries',{},'has_known',{},...
%         'loc_lines',{},'loc_cols',{},'loc_idx',{});
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Adrien Leygue (Adrien.Leygue@ec-nantes.fr)
%

result = WEAK([]);

    catalog_value = [];
    
    idx_res = 0;
    idx_value =0;
    for tt = 1:numel(WEAK)
        if WEAK(tt).has_known
            Nknown = numel(WEAK(tt).known_entries);
            fields_pos = old_catalog_value(WEAK(tt).known_entries,1);
            Nterms = [ALL_FIELDS(fields_pos).Nmodes];
            for i=1:prod(Nterms)
                idx_res = idx_res+1;
                result(idx_res) = WEAK(tt);
                tmp_catalog = old_catalog_value(WEAK(tt).known_entries,:);
                tmp_catalog(:,3) = my_ind2sub(Nterms,i);
                catalog_value(idx_value+(1:Nknown),:) = tmp_catalog; %#ok<AGROW>
                result(idx_res).known_entries = idx_value+(1:Nknown);
                idx_value = idx_value+Nknown;
            end
        else
            idx_res = idx_res+1;
            result(idx_res) = WEAK(tt);
        end
    end
    [catalog_value,~,scatter_value] = unique(catalog_value,'rows');
    RHS_mask = false(1,numel(result));
    for tt = 1:numel(result)
        RHS_mask(tt) = isempty(result(tt).unknown_entry);
        result(tt).known_entries = scatter_value(result(tt).known_entries);
    end
    RHS = result(RHS_mask);
    LHS = result(~RHS_mask);
end

function result = my_ind2sub(siz,idx)
result = zeros(numel(idx),numel(siz));
    k = [1 cumprod(siz(1:end-1))];
    for i = numel(siz):-1:1,
        vi = rem(idx-1, k(i)) + 1;
        vj = (idx - vi)/k(i) + 1;
        result(:,i) = vj;
        idx = vi;
    end
end