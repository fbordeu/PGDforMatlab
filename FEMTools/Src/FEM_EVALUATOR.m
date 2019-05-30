function [result, weights] = FEM_EVALUATOR(expr,DOM_NAME,COORDS,FIELDS,INTEG_PARAM,varargin)
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Adrien Leygue (Adrien.Leygue@ec-nantes.fr)
%
separated = false;
Nmodes = [FIELDS.Nmodes];
INTEG_VALUES = cell(numel(FIELDS),max([FIELDS.Ncomp])+1);
for i=1:numel(varargin)
    if ~ischar(varargin{i}), continue; end;
    switch upper(varargin{i})
        case 'SEPARATED'
            separated = varargin{i+1};
        case 'INTEG_VALUES'
            [found,pos] = ismember(varargin{i+1}{1},{FIELDS.name});
            assert(found,'provided fields not found among the known fields');
            assert(isempty(INTEG_VALUES{pos,1}),'duplicate provided field');
%            INTEG_VALUES{pos,numel(varargin{i+1})} = [];
            INTEG_VALUES(pos,1:(FIELDS(pos).Ncomp+1)) = varargin{i+1};
            Nmodes(pos) = size(varargin{i+1}{2},2);
            continue;
        otherwise
            error(['Unrecognized option:' varargin{i}]);
    end
end

%gather all fields and check if different fields (handles) have
%different names
ALL_FIELDS = unique([COORDS;FIELDS(:)]);
assert(all([ALL_FIELDS.type]~=2),'Value fields not supported');
assert(numel(unique({ALL_FIELDS.name}))==numel(ALL_FIELDS),'FEM_ASSEMBLE: Different fields share the same name');
ALL_FIELDS.set_current_domain(DOM_NAME);
COORDS.precompute_interpolants(DOM_NAME,INTEG_PARAM,ALL_FIELDS);



weights = COORDS.CURRENT_INTERP.PCIntegrationWeights{(COORDS.CURRENT_INTERP.PCIntegrationParam==INTEG_PARAM)};

if separated
    result = zeros(numel(weights),0);
else
    result = zeros(size(weights));
end
expr = process_expr(expr,FIELDS);
for t = 1:numel(expr)
    term = expr(t);
    if separated
        tmp = term.factor*ones([numel(weights) Nmodes(term.known_id)]);
    else
        tmp = term.factor*ones(size(weights));
    end
    for k = 1:numel(term.known_id)
        id = term.known_id(k);
        if (~isempty(INTEG_VALUES{id,term.known_comp(k)})) && (term.known_diff_symbol(k)==1)
            values = INTEG_VALUES{id,term.known_comp(k)+1};
        else
            INTERP = FIELDS(id).CURRENT_INTERP;
            values = (INTERP.PCMatrix{find(INTEG_PARAM==INTERP.PCIntegrationParam),term.known_diff_symbol(k)}*FIELDS(id).get_values('ALL',term.known_comp(k),1:Nmodes(id)));
        end
        if separated
            sz = [numel(weights) ones(1,numel(term.known_id))];
            sz(k+1) = Nmodes(id);
            tmp = bsxfun(@times,tmp,reshape(values,sz));
        else
            tmp = tmp.*values;
        end
    end
    if separated
        result = [result reshape(tmp,[numel(weights) prod(Nmodes(term.known_id))])];    
    else
    result = result+tmp;
    end
end
end

function expr = process_expr(expr,FIELDS)
for t = 1:numel(expr)
    kn_id = zeros(1,numel(expr(t).known_id));
    kn_diff = zeros(1,numel(expr(t).known_id));
    for k = 1:numel(expr(t).known_id)
        tmp = strcmp(expr(t).known_id{k},{FIELDS.name});
        msg = ['could not find ' expr(t).known_id{k} ' among the provided known fields'];
        assert(nnz(tmp)==1,msg);
        kn_id(k) = find(tmp);
        kn_diff(k) = ELEMENT.diff_id(expr(t).known_diff_symbol(:,k));
    end
    expr(t).known_id = kn_id;
    expr(t).known_diff_symbol = kn_diff;
end
end

