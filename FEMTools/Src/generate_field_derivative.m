function [DXYZDRST, varargout] =  generate_field_derivative(inField,fieldDerivative, compDerivative, spacename,GEOX,COORDS)
% inField is a vectoqr of  FEM_FIELD
% fieldDerivative is the space to apply the derivation
% compDerivative is the component of the derivations 
%
%                 1      2       fieldDerivative
%                 v      v
% for example u = f(x,y)*f(z)
%                   ^ ^    ^
%                   1 2    1     compDerivative
%
% if you want to  calculate du/dy :
%
% inField = u
% fieldDerivative = 1
% compDerivative = 2
% spacename = 'Y' ( to generate the name of the outputs fields)
% GEOX : a cell with the geo (to generate the correct interpolation)
% COORDS : a vector with the coordinates of the spaces 
%

%% Number of dims %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dims = numel(inField);


%% Generate the new fields %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%name of the new fields
name = [ 'D' inField(1).name(end:end) 'D' spacename ] ;
OOname = ['OO' name];

% interpolant of the new fields
for j = 1:dims
    interp_options = {};
    for i=1:numel(GEOX{j}.PHYSICAL_NAMES)
        interp_options = [interp_options {'sub' GEOX{j}.PHYSICAL_NAMES{i}}];  %#ok<AGROW>
    end
    %% two options with element with inferior order or with discontinues elements
    %%INTERP_CONST{j} = FEM_INTERP(GEOX{j},'ALL',GetDerivativeElement(inField(j)),interp_options);%#ok<AGROW>
    if fieldDerivative == j
        interp_options = [interp_options {'Discontinuous' }];  %#ok<AGROW> 
        INTERP_CONST{j} = FEM_INTERP(GEOX{j},'ALL',inField(j).INTERP(1).ELEM.name,interp_options);%#ok<AGROW>
    end
end


for i = 1: dims
    if i == fieldDerivative 
        DXYZDRST(i) = FEM_FIELD(name,INTERP_CONST{i});%#ok<AGROW>
        OODXYZDRST(i) = FEM_FIELD(OOname,INTERP_CONST{i});%#ok<AGROW>
    else
        DXYZDRST(i) = FEM_FIELD(name,inField(i).INTERP);%#ok<AGROW>
        OODXYZDRST(i) = FEM_FIELD(OOname,inField(i).INTERP);%#ok<AGROW>
    end
end

% first we copy all the terms, 
% make the derivation of the field
% clean the non relevant terms
% calculate the inverse

% copy all the terms
Ncomp = COORDS(fieldDerivative).Ncomp;
Fname = inField(i).name;
for i = 1 : dims
   if i == fieldDerivative
       coords = [];
       for c =1:Ncomp
           coords = [coords SYM_COORD(['coord' num2str(c)]) ];%#ok<AGROW>
       end
       SYM_M = SYM_FIELD(Fname,coords,'known',inField(i).Ncomp);% mapping
       SYM_D = SYM_FIELD(name,coords,'unknown',inField(i).Ncomp);% mapping derivative
       
       weak = SYM_D.test_fct'*(SYM_D-diff(SYM_M,coords(compDerivative)));
       [LHS, RHS] =  weak.extract_LHS_RHS();
       
       [A,b] = FEM_ASSEMBLE_PGD('ALL',COORDS(i),DXYZDRST(i),DXYZDRST(i),inField(i) ,2,'LHS',LHS.export,'RHS',RHS.export,'OLD_SCHOOL');
      
       X= A{1}\b;
       DXYZDRST(i).set_all_values(X)
   else
       DXYZDRST(i).set_all_values(inField(i).get_all_values())
   end
end

% Now we try to eliminate terms that share all the same terms in all dimension but one 
%D = cellfun(@(arg)(arg.VALUES),DXYZDRST ,'UniformOutput',false)';
D = {DXYZDRST.VALUES}';

nterms = size(D{1},2);
termsToKeep= true(1,nterms );

for i = 1:nterms -1
    for j = i+1:nterms 
        if ~termsToKeep(j) 
            continue
        end
        equal = true;
        dd = 0;
        for d = 1:dims
            if ~isequal(D{d}(:,i),D{d}(:,j))
                if equal == false % we accept only one dimension different
                    dd = 0;
                    break;
                end
                equal = false;
                dd = d;
            end
        end
        if dd > 0
            % addition of the terms 
            D{dd}(:,i) = D{dd}(:,i) + D{dd}(:,j);
            
            termsToKeep(j) = false;
        end
    end
end

% filter the terms to keep and put it back to the femfields
for d = 1:dims
    D{d} = D{d}(:,termsToKeep);
    DXYZDRST(d).set_all_values(D{d});
end

if nargout == 1 
    return 
end

% we calculate the inverse
sol = quotient(1,D,'max_added_modes',100*size(D{1},2));
sol = recompact(sol);

% clean empty modes

termsToKeep = sum(sol{1}'.^2,2)./norm(sol{1}(:,1)) > 1e-15;

for d = 1:dims
    sol{d} = sol{d}(:,termsToKeep);
    OODXYZDRST(d).set_all_values(sol{d});
end

varargout{1} =  OODXYZDRST;

end

function res= GetDerivativeElement(field)

name = field.INTERP(1).ELEM.name;

switch name
    case 'HEX_8' 
        res = 'HEX_1';
    case 'QUAD_4' 
        res = 'QUAD_1';
    case 'TRI_9' 
        res = 'TRI_3';
    case 'TRI_3' 
        res = 'TRI_1';
    case 'EDGE_2' 
        res = 'EDGE_1';
    otherwise
        error('Unexpected element type.');
end


end
%% 



