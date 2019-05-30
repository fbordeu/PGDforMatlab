function varargout = FEM_ASSEMBLE_PC(DOM_NAME,COORDS,U_FIELDS,T_FIELDS,K_FIELDS,INTEG_PARAM,varargin)
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Adrien Leygue (Adrien.Leygue@ec-nantes.fr)
%
LHS = [];
RHS = [];
if ~isempty(K_FIELDS);
INTEG_VALUES = cell(numel(K_FIELDS),max([K_FIELDS.Ncomp])+1);
end

for i=1:numel(varargin)
    if ~ischar(varargin{i}), continue; end;
    switch upper(varargin{i})
        case 'LHS'
            LHS = [LHS varargin{i+1}]; %#ok<AGROW>
            continue;
        case 'RHS'
            RHS = [RHS varargin{i+1}]; %#ok<AGROW>
            continue;
        case 'INTEG_VALUES'
            assert(~isempty(K_FIELDS),'INTEG_VALUES needs known fields');
            [found,pos] = ismember(varargin{i+1}{1},{K_FIELDS.name});
            assert(found,'provided fields not found among the known fields');
            assert(isempty(INTEG_VALUES{pos,1}),'duplicate provided field');
            %INTEG_VALUES{pos,numel(varargin{i+1})} = [];
            %INTEG_VALUES(pos,:) = varargin{i+1}; 
            INTEG_VALUES(pos,1:(K_FIELDS(pos).Ncomp+1)) = varargin{i+1};
            continue;
        otherwise
            error(['Unrecognized option:' varargin{i}]);
    end
end


%gather all fields and check if different fields (handles) have
%different names
ALL_FIELDS = unique([COORDS; T_FIELDS(:);U_FIELDS(:);K_FIELDS(:)]);
assert(all([ALL_FIELDS.type]~=2),'Value fields not supported');
assert(numel(unique({ALL_FIELDS.name}))==numel(ALL_FIELDS),'FEM_ASSEMBLE: Different fields share the same name');
ALL_FIELDS.set_current_domain(DOM_NAME);
COORDS.precompute_interpolants(DOM_NAME,INTEG_PARAM,ALL_FIELDS);

T_offsets = cumsum([T_FIELDS.Nnodes].*[T_FIELDS.Ncomp]);
T_offsets = T_offsets(:);
Nrows = T_offsets(end);
T_offsets = [0;T_offsets(1:end-1)];
T_comp_offsets = [T_FIELDS.Nnodes];

if ~isempty(U_FIELDS)
    U_offsets = cumsum([U_FIELDS.Nnodes].*[U_FIELDS.Ncomp]);
    U_offsets = U_offsets(:);
    Ncols = U_offsets(end);
    U_offsets = [0;U_offsets(1:end-1)];
    U_comp_offsets = [U_FIELDS.Nnodes];
end

weights = COORDS.CURRENT_INTERP.PCIntegrationWeights{(COORDS.CURRENT_INTERP.PCIntegrationParam==INTEG_PARAM)};

has_rhs = ~isempty(RHS);
has_lhs = ~isempty(LHS);
if has_lhs
    LHS = LHS(find([LHS.factor]));
    all_i = [];
    all_j = [];
    all_v = [];
    LHS = process_weak(LHS,U_FIELDS,T_FIELDS,K_FIELDS);
    for t = 1:numel(LHS)
        term = LHS(t);
        tmp = term.factor*weights;
        U_INTERP = U_FIELDS(term.unknown_id).CURRENT_INTERP;
        T_INTERP = T_FIELDS(term.test_id).CURRENT_INTERP;
        for k = 1:numel(term.known_id)
            id = term.known_id(k);
            if (~isempty(INTEG_VALUES{id,term.known_comp(k)})) && (term.known_diff_symbol(k)==1)
                tmp = tmp.* INTEG_VALUES{id,term.known_comp(k)+1}(:);
            else
                K_INTERP = K_FIELDS(id).CURRENT_INTERP;
                tmp = tmp.* (K_INTERP.PCMatrix{find(INTEG_PARAM==K_INTERP.PCIntegrationParam),term.known_diff_symbol(k)}*K_FIELDS(id).get_values('ALL',term.known_comp(k))); %#ok<*FNDSB>
            end
        end
        %tmpA = bsxfun(@times,T_INTERP.PCMatrix{find(INTEG_PARAM==T_INTERP.PCIntegrationParam),term.test_diff_symbol},tmp)'*U_INTERP.PCMatrix{find(INTEG_PARAM==U_INTERP.PCIntegrationParam),term.unknown_diff_symbol};
        tmpA = T_INTERP.PCMatrix{find(INTEG_PARAM==T_INTERP.PCIntegrationParam),term.test_diff_symbol}'*spdiags(tmp,0,numel(tmp),numel(tmp))*U_INTERP.PCMatrix{find(INTEG_PARAM==U_INTERP.PCIntegrationParam),term.unknown_diff_symbol};

        [tmp_i,tmp_j,tmp_v] = find(tmpA);
        tmp_i = T_offsets(term.test_id)+ T_comp_offsets(term.test_id)*(term.test_comp-1)+tmp_i;
        tmp_j = U_offsets(term.unknown_id)+ U_comp_offsets(term.unknown_id)*(term.unknown_comp-1)+tmp_j;
        
        all_i = [all_i;tmp_i]; %#ok<AGROW>
        all_j = [all_j;tmp_j]; %#ok<AGROW>
        all_v = [all_v;tmp_v]; %#ok<AGROW>
    end
    A = sparse(all_i,all_j,all_v,Nrows,Ncols);
end
if has_rhs
    RHS = RHS(find([RHS.factor]));
    b = zeros(Nrows,1);
    RHS = process_weak(RHS,U_FIELDS,T_FIELDS,K_FIELDS);
    for t = 1:numel(RHS)
        term = RHS(t);
        T_INTERP = T_FIELDS(term.test_id).CURRENT_INTERP;
        tmp = term.factor*weights;
        for k = 1:numel(term.known_id)
            id = term.known_id(k);
            if (~isempty(INTEG_VALUES{id,term.known_comp(k)})) && (term.known_diff_symbol(k)==1)
                tmp = tmp.* INTEG_VALUES{id,term.known_comp(k)+1}(:);
            else
                
                K_INTERP = K_FIELDS(id).CURRENT_INTERP;
                tmp = tmp.* (K_INTERP.PCMatrix{find(INTEG_PARAM==K_INTERP.PCIntegrationParam),term.known_diff_symbol(k)}*K_FIELDS(id).get_values('ALL',term.known_comp(k)));
            end
        end
        tmpb = tmp'* T_INTERP.PCMatrix{find(INTEG_PARAM==T_INTERP.PCIntegrationParam),term.test_diff_symbol};
        b_offset = T_offsets(term.test_id)+ T_comp_offsets(term.test_id)*(term.test_comp-1);
        b((1:numel(tmpb))+b_offset) = b((1:numel(tmpb))+b_offset) + tmpb';
    end
end

varargout = cell(1,2);
varargout{1} = sparse(Nrows,Ncols);
varargout{2} = sparse(Nrows,1);
if has_lhs
    varargout{1} = A;
end
if has_rhs
    varargout{2} = b;
end
if nargout==1
    varargout = varargout([has_lhs has_rhs]);
end
end

function weak = process_weak(weak,U_FIELDS,T_FIELDS,K_FIELDS)
for t = 1:numel(weak)
    tmp = strcmp(weak(t).test_id{1},{T_FIELDS.name});
    msg = ['could not find ' weak(t).test_id{1} ' among the provided test fields'];
    assert(nnz(tmp)==1,msg);
    weak(t).test_id = find(tmp);
    weak(t).test_diff_symbol = ELEMENT.diff_id(weak(t).test_diff_symbol);
    
    if ~isempty(weak(t).unknown_id)
        tmp = strcmp(weak(t).unknown_id{1},{U_FIELDS.name});
        msg = ['could not find ' weak(t).unknown_id{1} ' among the provided unknown fields'];
        assert(nnz(tmp)==1,msg);
        weak(t).unknown_id = find(tmp);
        weak(t).unknown_diff_symbol = ELEMENT.diff_id(weak(t).unknown_diff_symbol);
    end
    
    kn_id = zeros(1,numel(weak(t).known_id));
    kn_diff = zeros(1,numel(weak(t).known_id));
    for k = 1:numel(weak(t).known_id)
        tmp = strcmp(weak(t).known_id{k},{K_FIELDS.name});
        msg = ['could not find ' weak(t).known_id{k} ' among the provided known fields'];
        assert(nnz(tmp)==1,msg);
        kn_id(k) = find(tmp);
        kn_diff(k) = ELEMENT.diff_id(weak(t).known_diff_symbol(:,k));
    end
    weak(t).known_id = kn_id;
    weak(t).known_diff_symbol = kn_diff;
end
end

