function [AII,AIV,AJJ,AJV,AW,Aoffsets,BII,BIV,BW,Boffsets,Nrows,Ncols] = PGD_DETAILLED_ASSEMBLE(DOM_NAME,COORDS,U_FIELDS,T_FIELDS,K_FIELDS,INTEG_PARAM,varargin)
LHS = [];
RHS = [];
if ~isempty(K_FIELDS)
INTEG_VALUES = cell(numel(K_FIELDS),max([K_FIELDS.Ncomp]));
K_Nmodes = [K_FIELDS.Nmodes];
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
            [found,pos] = ismember(varargin{i+1}{1},{K_FIELDS.name});
            assert(found,'provided fields not found among the known fields');
            assert(isempty(INTEG_VALUES{pos,1}),'duplicate provided field');
            INTEG_VALUES{pos,numel(varargin{i+1})} = [];
            INTEG_VALUES(pos,:) = varargin{i+1};
            K_Nmodes(pos) = size(varargin{i+1}{2},2);
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

Ncols = 1;
if ~isempty(U_FIELDS)
    U_offsets = cumsum([U_FIELDS.Nnodes].*[U_FIELDS.Ncomp]);
    U_offsets = U_offsets(:);
    Ncols = U_offsets(end);
    U_offsets = [0;U_offsets(1:end-1)];
    U_comp_offsets = [U_FIELDS.Nnodes];
end

weights = COORDS.CURRENT_INTERP.PCIntegrationWeights{(COORDS.CURRENT_INTERP.PCIntegrationParam==INTEG_PARAM)};
Nweights = numel(weights);

if(~isempty(LHS))
    LHS = LHS(find([LHS.factor]));
end
if (~isempty(RHS))
    RHS = RHS(find([RHS.factor]));
end

has_rhs = ~isempty(RHS);
has_lhs = ~isempty(LHS);

if has_lhs
    
    LHS = process_weak(LHS,U_FIELDS,T_FIELDS,K_FIELDS);
    gen_terms= ones(1,numel(LHS));
    if ~isempty(K_FIELDS)
    for t = 1:numel(LHS)
        gen_terms(t) = prod(K_Nmodes(LHS(t).known_id));
    end
    end
    AII = cell(1,sum(gen_terms));
    AIV= cell(1,sum(gen_terms));
    AJJ= cell(1,sum(gen_terms));
    AJV= cell(1,sum(gen_terms));
    AW= cell(1,sum(gen_terms));
    Aoffsets = cell(1,sum(gen_terms));
    
    idx_AA = 1;
    for t = 1:numel(LHS)
        term = LHS(t);
        tmp = term.factor*weights;
        U_INTERP = U_FIELDS(term.unknown_id).CURRENT_INTERP;
        T_INTERP = T_FIELDS(term.test_id).CURRENT_INTERP;
        for k = 1:numel(term.known_id)
            id = term.known_id(k);
            if (~isempty(INTEG_VALUES{id,term.known_comp(k)})) && (term.known_diff_symbol(k)==1)
                values = INTEG_VALUES{id,term.known_comp(k)+1};
            else
            K_INTERP = K_FIELDS(id).CURRENT_INTERP;
            values = K_INTERP.PCMatrix{find(INTEG_PARAM==K_INTERP.PCIntegrationParam),term.known_diff_symbol(k)}*K_FIELDS(id).get_values('ALL',term.known_comp(k),1:K_FIELDS(id).Nmodes);
            end
            tmp = reshape(bsxfun(@times,reshape(tmp,[Nweights 1 size(tmp,2)]),values),Nweights,size(tmp,2)*size(values,2));
        end
        AT = T_INTERP.PCMatrix{find(INTEG_PARAM==T_INTERP.PCIntegrationParam),term.test_diff_symbol};
        [aii(:,1),aii(:,2),aiv] = find(AT);
        AU = U_INTERP.PCMatrix{find(INTEG_PARAM==U_INTERP.PCIntegrationParam),term.unknown_diff_symbol};
        [ajj(:,1),ajj(:,2),ajv] = find(AU);
        for i=1:gen_terms(t)
            AII{idx_AA} =  aii;
            AIV{idx_AA} = aiv;
            AJJ{idx_AA} = ajj;
            AJV{idx_AA} = ajv;
            AW{idx_AA}  = tmp(:,i);
            Aoffsets{idx_AA} = [T_offsets(term.test_id)+T_comp_offsets(term.test_id)*(term.test_comp-1)...
                                U_offsets(term.unknown_id)+U_comp_offsets(term.unknown_id)*(term.unknown_comp-1) ];
            idx_AA = idx_AA+1;
        end
    end
    
else
    AII = cell(1,0);
    AIV= cell(1,0);
    AJJ= cell(1,0);
    AJV= cell(1,0);
    AW= cell(1,0);
    Aoffsets = cell(1,0); 
end

if has_rhs
    RHS = process_weak(RHS,U_FIELDS,T_FIELDS,K_FIELDS);
    gen_terms= ones(1,numel(RHS));
    if ~isempty(K_FIELDS)
    for t = 1:numel(RHS)
        gen_terms(t) = prod(K_Nmodes(RHS(t).known_id));
    end
    end
    BII = cell(1,numel(RHS));
    BIV= cell(1,numel(RHS));
    BW= cell(1,numel(RHS));
    Boffsets = cell(1,numel(RHS));

    for t = 1:numel(RHS)
        term = RHS(t);
        T_INTERP = T_FIELDS(term.test_id).CURRENT_INTERP;
        tmp = term.factor*weights;
        for k = 1:numel(term.known_id)
            id = term.known_id(k);
            if (~isempty(INTEG_VALUES{id,term.known_comp(k)})) && (term.known_diff_symbol(k)==1)
                values = INTEG_VALUES{id,term.known_comp(k)+1};
            else
            K_INTERP = K_FIELDS(id).CURRENT_INTERP;
            values = K_INTERP.PCMatrix{find(INTEG_PARAM==K_INTERP.PCIntegrationParam),term.known_diff_symbol(k)}*K_FIELDS(id).get_values('ALL',term.known_comp(k),1:K_FIELDS(id).Nmodes);
            end
            tmp = reshape(bsxfun(@times,reshape(tmp,[Nweights 1 size(tmp,2)]),values),Nweights,size(tmp,2)*size(values,2));
        end
        [bii(:,1),bii(:,2),biv] = find(T_INTERP.PCMatrix{find(INTEG_PARAM==T_INTERP.PCIntegrationParam),term.test_diff_symbol}); %#ok<*FNDSB>
        BII{t} =bii;
        BIV{t} = biv;
        BW{t}  = tmp;
        Boffsets{t} = T_offsets(term.test_id)+ T_comp_offsets(term.test_id)*(term.test_comp-1);
    end
else
    BII = cell(1,0);
    BIV= cell(1,0);
    BW= cell(1,0);  
    Boffsets = cell(1,0);
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

