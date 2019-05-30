%SYMBOLIC PDE TOOLBOX BASE CLASS
%
%Abstract base class for the PDE symbolic toolbox
%All operations on symbolic entities go through this class which dispaches
%specialized operation to protected methods of the derivec classes.
%
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Adrien Leygue (Adrien.Leygue@ec-nantes.fr)
%

classdef (Abstract = true)  SYM_BASE < matlab.mixin.Heterogeneous
    properties
        %string identifying the CLASS
        name = '';
        %true if the class represents a real number with a determined value
        num = false;
    end
    
    
    
    methods (Sealed = true)
        
        %EXPORT & DISPLAY FUNCTIONS
        function display(arg)
            if ~isempty(inputname(1))
                display([' ' inputname(1) '=']);
            end
            display(sym2str(arg));
        end
        
        function disp(arg)
            if ~isempty(inputname(1))
                display([' ' inputname(1) '=']);
            end
            display(sym2str(arg));
        end
        
        
        function result = sym2str(arg)
            %Converts a symbolic expression to a string
            %For non-scalar arguments, a cell array with the string of each
            %element is returned
            if isscalar(arg)
                result = sym2str_private(arg);
            else
                sz = size(arg);
                result = cell(sz);
                for i=1:numel(arg)
                    result{i} = sym2str_private(arg(i));
                end
            end
        end
        
        function result = sym2double(arg)
            result = zeros(size(arg));
            for i=1:numel(arg)
                result(i) = sym2double_private(arg(i));
            end
        end
        
        function result = display_tex(arg)
            %Converts a symbolic expression to LaTeX string.
            %For non-scalar arguments, a cell array with the string of each
            %element is returned
            if isscalar(arg)
                result = display_tex_private(arg);
            else
                result = cell(size(arg));
                for i=1:numel(arg)
                    result{i} = display_tex_private(arg(i));
                end
            end
        end
        %function result = FEM_export(arg)
        function result = export(arg)
            
            Nterms = numel(arg);
            result = struct('factor',{},...
                'test_id',{},'test_comp',{},'test_diff_symbol',{},...
                'unknown_id',{},'unknown_comp',{},'unknown_diff_symbol',{},...
                'known_id',{},'known_comp',{},'known_diff_symbol',{},'known_depend_symbol',{});
            for t = Nterms:-1:1
                F = get_factors_private(arg(t));
                result(t).factor = 1;
                for i=1:numel(F)
                    if F(i).num
                        result(t).factor = result(t).factor*F(i).value;
                    elseif isa(F(i),'SYM_FIELD')
                        switch F(i).tag
                            case 'test'
                                result(t).test_id = {F(i).id};
                                result(t).test_comp = F(i).comp;
                                result(t).test_diff_symbol = F(i).diff_symbol(:);
                            case 'unknown'
                                result(t).unknown_id = {F(i).id};
                                result(t).unknown_comp = F(i).comp;
                                result(t).unknown_diff_symbol = F(i).diff_symbol(:);
                            case 'known'
                                result(t).known_id = [ result(t).known_id {F(i).id}];
                                result(t).known_comp = [ result(t).known_comp F(i).comp];
                                result(t).known_diff_symbol = [ result(t).known_diff_symbol F(i).diff_symbol(:)];
                                result(t).known_depend_symbol = [ result(t).known_depend_symbol F(i).depend_symbol(:)];
                                
                            otherwise
                                error('unexpected tag');
                        end
                    else
                        error('SYM_BASE: FEM_export only SYM_FIELDS and Numerical types should be present ar this point');
                    end
                    
                end
            end
        end
        
        
        
        function result = remove_products_of_unknown_fields(arg)
            result = SYM_ZERO(size(arg));
            for i=1:numel(arg)
                result(i) = remove_products_of_unknown_fields_private(arg(i));
            end
        end
        
        function [LHS,RHS] = extract_LHS_RHS(arg)
            %Splits a PDE in weighted residual form into the left and
            %right hand side
            %LHS contains the different terms that appear in the Left Hand Side
            %RHS contains the different terms that appear in the Right Hand Side.
            
            %Extracts all the terms of the expression
            if isa(arg,'SYM_SUM')
                T = arg.terms;
            else
                T = arg;
            end
            %basic check for a test function
            assert(all(contains_tagged_field(T,'test')),'SYM_BASE: extract_LHS_RHS: All terms in the expression should contain a test function');
            mask = contains_tagged_field(T,'unknown');
            LHS = SYM_ZERO; % default value
            RHS = SYM_ZERO; % default value
            if any(mask)
                LHS = T(mask);
            end
            if any(~mask)
                RHS = -T(~mask);
            end
        end
        
        function varargout = export_expression(arg,opt)
            %Splits an expression into its different terms
            %The expression should not contain any test function
            %all fields are treated as known fields
            
            %if the expression is not scalar returns one output per element
            varargout = cell(1,numel(arg));
            if ~isscalar(arg)
                for i=1:numel(arg)
                    varargout{i} = export_expression(arg(i));
                end
                if nargin>1
                    if strcmp(opt,'nocell')
                        if numel(varargout)>1
                            assert(all(cellfun(@isscalar,varargout)),'nocell option not applicable, sum in the expression');
                            tmp = cell2mat(varargout);
                            varargout = {reshape(tmp,size(arg))};
                        end
                    else
                        error('unknown option')
                    end
                end
                
                return;
            end
            %Extracts all the terms of the expression
            if isa(arg,'SYM_SUM')
                T = arg.terms;
            else
                T = arg;
            end
            Nterms = numel(T);
            %basic check for a test function
            assert(~any(contains_tagged_field(T,'test')),'SYM_BASE: extract_expression: test functions found in the expression');
            result = struct('factor',{},...
                'known_id',{},'known_comp',{},'known_diff_symbol',{},'known_depend_symbol',{});
            for t = Nterms:-1:1
                F = get_factors_private(T(t));
                result(t).factor = 1;
                for i=1:numel(F)
                    if F(i).num
                        result(t).factor = result(t).factor*F(i).value;
                    elseif isa(F(i),'SYM_FIELD')
                        result(t).known_id = [ result(t).known_id {F(i).id}];
                        result(t).known_comp = [ result(t).known_comp F(i).comp];
                        result(t).known_diff_symbol = [ result(t).known_diff_symbol F(i).diff_symbol(:)];
                        result(t).known_depend_symbol = [ result(t).known_depend_symbol F(i).depend_symbol(:)];
                    else
                        error('SYM_BASE: FEM_export only SYM_FIELDS and Numerical types should be present ar this point');
                    end
                end
            end
            varargout{1} = result;
        end
        
        function varargout = separate(arg,pattern)
            if ~isscalar(arg)
                for i=1:numel(arg)
                    varargout{i} = separate(arg(i),pattern);
                end
                return;
            end
            %              if isa(arg,'SYM_SUM')
            %                 T = arg.terms;
            %             else
            %                 T = arg;
            %             end
            %             Nterms = numel(T);
            %             for i=1:Nterms
            %
            result = separate_private(arg,pattern);
            result(2,:) = result(1,:).*result(2,:);
            result = result(2:end,:);
            varargout = {result};
        end
        
        function [LHS,RHS,LHS_groups,RHS_groups] = extract_separated_LHS_RHS(arg,pattern)
            %Splits a PDE in weighted residual form into the left and
            %right hand side, and performs a separation of variable
            %according to the provided pattern (partition of the
            %coordinates)
            %pattern has the following form: { [coord1 coord2], coord3, [coord4]}
            %LHS contains the different terms that appear in the Left Hand Side
            %The columns of LHS contains the different terms that appear in the Left Hand Side
            %The Columns of RHS contains the different terms that appear in the Right Hand Side.
            %for each column of the outputs, the different lines correspond
            %to the separated contribution along each separated space.
            %NB: The first line contains the scalar factors.
            
            %split before separation of variables
            [LHS_tmp,RHS_tmp] = extract_LHS_RHS(arg);
            for i=1:numel(LHS_tmp)
                %call to class specific implementation
                LHS(:,i) = separate_private(LHS_tmp(i),pattern); %#ok<AGROW>
            end
            for i=1:numel(RHS_tmp)
                %call to class specific implementation
                RHS(:,i) = separate_private(RHS_tmp(i),pattern); %#ok<AGROW>
            end
            LHS(2,:) = LHS(2,:).*LHS(1,:);
            RHS(2,:) = RHS(2,:).*RHS(1,:);
            LHS = LHS(2:end,:);
            RHS = RHS(2:end,:);
            LHS_groups = 1:size(LHS,2);
            RHS_groups = 1:size(RHS,2);
            
        end
        
        %ADVANCED MANIPULATION (DIFF, SUBS,...)
        function result = subs(arg,old,new)
            %Substitute a symbolic expression with another expression
            %RESULT  = subs(ARG,OLD,NEW)
            %in ARG, all occurences of OLD are replaced by NEW
            %if OLD and NEW are arrays, replaces OLD(1) by ARG(1), then
            %uses this result to replaces OLD(2), by ARG(2),....
            %if NEW is numeric, it is converted to a symbolic object first
            
            new = conv2SYM(new);
            assert(numel(old)==numel(new),'SYM_BASE: subs needs old and new to be the same size')
            result = arg;
            for k=1:numel(old)
                for i=1:numel(result)
                    %call to class specific implementation
                    result(i) = subs_private(result(i),old(k),new(k));
                end
            end
        end
        
        function result = test_fct(arg)
            %Changes an unknown function to a test function
            %works as it should for more complex expressions
            for i=1:numel(arg)
                %call to class specific implementation
                result(i) = test_fct_private(arg(i)); %#ok<AGROW>
            end
            result = reshape(result,size(arg));
        end
        
        function result = diff(arg,var)
            %computes the derivative of arg with respect to the coordinate var
            for i=numel(arg):-1:1
                result(i) = diff_private(arg(i),var);
            end
            result = reshape(result,size(arg));
        end
        
        function result = divergence(arg,vars)
            %computes the divergence of arg with respect to the coordinates vars (array of SYM_COORD)
            %if arg is  multidimensional array, the function operates along the
            %last dimension of dimension numel(vars)
            % div(a_ij) = D a_ij / D x_j
            
            sz = size(arg);
            my_dim = find(sz==numel(vars),1,'last');
            assert(~isempty(my_dim),'SYM_BASE: divergence: incompatible dimensions')
            
            result = SYM_ZERO(sz);
            for i=1:prod(sz)
                my_sub = my_ind2sub(sz,i);
                result(i) = diff(arg(i),vars(my_sub(my_dim)));
            end
            result = sum(result,my_dim);
            result = squeeze(result);
        end
        
        function result  = gradient(arg,vars)
            %computes the gradient of a with respect to the coordinates
            %vars (array of SYM_COORD). The function will add a new
            %dimension of size numel(vars) and then squeeze any
            %singleton dimension
            %size(result) = squeeze([size(arg) numel(vars)])
            %gradient(U_i) = D U_i / D x_j
            result(numel(arg)*numel(vars)) = SYM_ZERO;
            idx = 1;
            for i=1:numel(vars)
                for j=1:numel(arg)
                    result(idx) = diff(arg(j),vars(i));
                    idx = idx+1;
                end
            end
            tmp = [size(arg) numel(vars)];
            result = squeeze(reshape(result,tmp));
        end
        
        function result = curl(arg,vars)
            %computes the curl of the vector ARG with respect tothe
            %variables VARS
            assert(numel(arg)==3,'SYM_BASE: curl applies to vectors of length 3');
            assert(numel(vars)==3,'SYM_BASE: curl needs 3 variables do differentiate');
            result = SYM_ZERO(size(arg));
            result(1) = diff(arg(3),vars(2)) - diff(arg(2),vars(3));
            result(2) = diff(arg(1),vars(3)) - diff(arg(3),vars(1));
            result(3) = diff(arg(2),vars(1)) - diff(arg(1),vars(2));
        end
        
        function result = sum(arg,opt)
            %overloads the sum function
            if nargin==1
                opt = 1;
            end
            sz = size(arg);
            if isscalar(sz)
                result = arg;
                return;
            end
            new_sz = sz;
            new_sz(opt) = 1;
            result = SYM_ZERO(new_sz);
            for i=1:numel(arg)
                pos = cell(size(sz));
                [pos{:}] = ind2sub(sz,i);
                pos{opt} = 1;
                %could be improved by a call to atomic times and a final
                %collect_terms on  result
                result(pos{:}) = atomic_plus(result(pos{:}),arg(i));
            end
            for i=1:numel(result)
                if isa(result(i),'SYM_SUM')
                    result(i) = collect_scalars(result(i));
                end
            end
        end
        
        function result = prod(arg,opt)
            %overloads the prod function
            if nargin==1
                opt = 1;
            end
            sz = size(arg);
            if isscalar(sz)
                result = arg;
                return;
            end
            new_sz = sz;
            new_sz(opt) = 1;
            result = SYM_ONE(new_sz);
            for i=1:numel(arg)
                pos = cell(size(sz));
                [pos{:}] = ind2sub(sz,i);
                pos{opt} = 1;
                result(pos{:}) = atomic_times(result(pos{:}),arg(i));
            end
        end
        
        function result = diag(obj,k)
            %overloads the diag function
            if nargin==1
                k=0;
            end
            mod_line = max(0,-k);
            mod_row = max(0,k);
            if isvector(obj)
                result = SYM_ZERO(numel(obj)+k);
                for i= 1:numel(obj)
                    result(i+mod_line,i+mod_row) = obj(i);
                end
            else
                sz = size(obj);
                assert(numel(sz)==2&(sz(1)==sz(2)),'SYM_BASE: diag needs a square matrix')
                result = SYM_ZERO([sz(1)-abs(k) 1]);
                for i=1:(sz(1)-abs(k))
                    result(i) = obj(i+mod_line,i+mod_row);
                end
            end
        end
        
        function result = trace(arg)
            %overloads trace
            result = sum(diag(arg));
        end
        
        function result = inv(arg)
            %overloads inv for scalar quantities through power
            assert(isscalar(arg),'SYM_BASE: inv not implemented for non scalar inputs');
            result = inv_private(arg);
        end
        
        function result = sin(arg)
            %overloads sin
            result = SYM_ZERO(size(arg));
            for i=1:numel(arg)
                result(i) = SYM_FUNCTION('sin',arg(i));
            end
        end
        
        function result = cos(arg)
            %overloads cos
            result = SYM_ZERO(size(arg));
            for i=1:numel(arg)
                result(i) = SYM_FUNCTION('cos',arg(i));
            end
        end
        
        function result = exp(arg)
            %overloads exp
            result = SYM_ZERO(size(arg));
            for i=1:numel(arg)
                result(i) = SYM_FUNCTION('exp',arg(i));
            end
        end
        
        %BASIC OPERATOR OVERLOAD
        
        function result = power(a,b)
            a = conv2SYM(a);
            b = conv2SYM(b);
            [a,b] = match_scalar(a,b);
            for i=numel(a):-1:1
                %call to non-vectorized function
                %result(i) = simplify_private(SYM_FUNCTION('power',a(i),b(i)));
                result(i) = (SYM_FUNCTION('power',a(i),b(i)));
                
            end
            result = reshape(result,size(a));
        end
        
        function result = sqrt(a)
            result = power(a,SYM_REAL(1/2));
        end
        
        function result = mpower(a,b)
            if(isscalar(b))
                result=power(a,b);
                return;
            end
            error('SYM_BASE: mpower not implemented for non scalar second argument');
        end
        
        function result = plus(a,b)
            a = conv2SYM(a);
            b = conv2SYM(b);
            [a,b] = match_scalar(a,b);
            for i=numel(a):-1:1
                %call to non-vectorized function
                result(i) = atomic_plus(a(i),b(i));
                if isa(result(i),'SYM_SUM')
                    result(i) = collect_scalars(result(i));
                end
            end
            result = reshape(result,size(a));
        end
        
        function result = times(a,b)
            a = conv2SYM(a);
            b = conv2SYM(b);
            [a,b] = match_scalar(a,b);
            for i=numel(a):-1:1
                %call to non-vectorized function
                result(i) = atomic_times(a(i),b(i));
            end
            result = reshape(result,size(a));
        end
        
        function result=mtimes(a,b)
            if(isscalar(a)||isscalar(b))
                result=times(a,b);
                return;
            end
            a = conv2SYM(a);
            b = conv2SYM(b);
            assert((size(a,2)==size(b,1)),'SYM_BASE: mtimes(a,b): a and b must have compatible sizes');
            m = size(a,1);
            n = size(b,2);
            o = size(a,2);
            for i=m:-1:1
                for j = n:-1:1
                    %call to non-vectorized function
                    result(i,j) = atomic_times(a(i,1),b(1,j));
                    for k=2:o
                        %call to non-vectorized function
                        result(i,j) = atomic_plus(result(i,j),atomic_times(a(i,k),b(k,j)));
                    end
                    if isa(result(i,j),'SYM_SUM')
                        result(i,j) = collect_scalars(result(i,j));
                    end
                end
            end
        end
        
        function result=mrdivide(a,b)
            if(isscalar(b))
                result=rdivide(a,b);
                return;
            end
            error('SYM_BASE: mrdivide not implemented for non scalar second argument');
        end
        
        function result = rdivide(a,b)
            a = conv2SYM(a);
            b = conv2SYM(b);
            [a,b] = match_scalar(a,b);
            for i=numel(a):-1:1
                %call to non-vectorized function
                result(i) = atomic_rdivide(a(i),b(i));
            end
            result = reshape(result,size(a));
        end
        
        function result = uminus(a)
            result = times(a,SYM_REAL(-1));
        end
        
        function result=minus(a,b)
            result = plus(a,uminus(b));
        end
        
        function result=uplus(a)
            result = a;
        end
        
        function result = eq(a,b)
            if isscalar(a)
                result = true(size(b));
                for i=1:numel(b)
                    %call to non-vectorized function
                    result(i) = atomic_eq(a,b(i));
                end
                return;
            end
            if isscalar(b)
                result = true(size(a));
                for i=1:numel(a)
                    %call to non-vectorized function
                    result(i) = atomic_eq(a(i),b);
                end
                return;
            end
            assert(all(size(a)==size(b)),'SYM_BASE:eq(a,b) if a and b are not scalars,they must have the same size')
            result = true(size(b));
            for i=1:numel(b)
                %call to non-vectorized function
                result(i) = atomic_eq(a(i),b(i));
            end
        end
        
        function result = simplify_products(arg,opt)
            if nargin==1
                opt='collect';
            end
            for i=numel(arg):-1:1
                result(i) = simplify_products_private(arg(i),opt);
            end
            result = reshape(result,size(arg));
        end
    end % SEALED METHODS
    
    methods ( Sealed, Access = protected)
        %INTERNAL UTILITY FUNCTIONS
        function tf = contains_tagged_field(arg,tag)
            %true if arg contains a SYM_FIELD object with tag TAG
            %false otherwise
            tf = false(size(arg));
            for i=1:numel(arg)
                %call to class specific implementation
                tf(i) = contains_tagged_field_private(arg(i),tag);
            end
        end
        
        %ATOMIC FUNCTIONS FOR OPERATORS
        function result = atomic_eq(a,b)
            %eq function for two scalar arguments
            if strcmp(class(a),class(b))
                %call to class specific implementation
                result = eq_private(a,b);
            else
                result = false;
            end
        end
        
        function result = atomic_plus(a,b)
            %times function for scalar arguments
            
            %Particular cases for which the output is one of the arguments
            if(isa(a,'SYM_ZERO'))
                result = b;
                return
            end
            if(isa(b,'SYM_ZERO'))
                result = a;
                return
            end
            %in case we have two numeric arguments
            if(a.num&&b.num)
                result = SYM_REAL(a.value+b.value);
                result = atomic_simplify_num(result);
            else
                %General case, the result is simplified afterwards
                result = SYM_SUM(a,b);
                %result = collect_scalars(result);
            end
        end
        
        function result = atomic_times(a,b)
            %times function for scalar arguments
            
            %Particular cases for which rhe result is one of the arguments
            if(isa(a,'SYM_ZERO'))
                result = a;
                return
            end
            if(isa(b,'SYM_ZERO'))
                result = b;
                return
            end
            if(isa(a,'SYM_ONE'))
                result = b;
                return
            end
            if(isa(b,'SYM_ONE'))
                result = a;
                return
            end
            
            %in case we have two numeric arguments
            if(a.num&&b.num)
                result = SYM_REAL(a.value*b.value);
                result = atomic_simplify_num(result);
                return;
            end
            
            %for a sum, one has to distribute
            if (isa(a,'SYM_SUM'))
                result = SYM_ZERO;
                for i=1:numel(a.terms)
                    result = atomic_plus(result,atomic_times(a.terms(i),b));
                end
                return;
            end
            if (isa(b,'SYM_SUM'))
                result = SYM_ZERO;
                for i=1:numel(b.terms)
                    result = atomic_plus(result,atomic_times(b.terms(i),a));
                end
                return;
            end
            %general case
            result = SYM_PROD(a,b);
            
        end
        
        function result = atomic_rdivide(a,b)
            % rdivide for scalar arguments
            if isa(b,'SYM_ZERO')
                error('Cannot divide by SYM_ZERO');
            end
            result = atomic_times(a,inv(b));
        end
        
        
        function [result,str,p] = sort_disp(arg)
            %sorts the elements of arg according to the order generated by
            %their sym2str string.
            %result is the sorted vector
            %str is the cell array of sorted strings
            %p is the sort permutation
            %str = sym2str(arg);
            str = cell(1,numel(arg));
            for i=1:numel(arg)
                str{i} = sym2str_private(arg(i));
            end
            [str,p] = sort(str);
            result = arg(p);
        end
        
        
        function [scalars,entities] = normalize(arg)
            %normalizes an array to an array of real numbers and an array of normalized
            %entities.
            for i=1:numel(arg)
                %call to the class specific implementation
                [scalars(i), entities(i)] = normalize_private(arg(i)); %#ok<AGROW>
            end
        end
    end
    
    methods ( Hidden = true, Sealed)
        function result = create_handle_function(arg,coords)
            str = ['@(' sym2str_private(coords(1))];
            for i=2:numel(coords)
                str = [str ',' sym2str(coords(i))]; %#ok<AGROW>
            end
            if isscalar(arg)
                str = [str ')' sym2str(arg)];
            else
                tmp = sym2str(arg);
                str = [str ') reshape(['  strjoin(tmp(:),' , ') '],' num2str(size(arg,1)) ',' num2str(size(arg,2)) ')'];
            end
            result = eval(str);
        end
        
        function result = atomic_simplify_num(arg)
            if arg.num
                switch arg.value
                    case 0
                        result = SYM_ZERO;
                    case 1
                        result = SYM_ONE;
                    otherwise
                        result = arg;
                end
            else
                result = arg;
            end
        end
        
    end
    
    
    
    
    %CLASS SPECIFIC FUNCTIONS WITH NO DEFAULT IMPLEMENTATION
    methods (Abstract = true, Access = protected)
        %THESE METHODS HAVE TO BE IMPLEMENTED IN EACH CASE
        result = display_tex_private(arg)
        result = sym2str_private(arg)
        result = eq_private(a,b)
        result = export_private(arg)
        result = separate_private(arg,pattern);
        
    end
    
    
    %DEFAULT IMPLEMENTATION
    methods (Access = protected)
        %THESE FUNCTIONS MIGHT BE OVERLOADED -->NOT SEALED
        function tf = contains_tagged_field_private(~,~)
            %checks if arg contains a SYM_FIELD tagged tag.
            %default implementation.
            %this function is overloaded for specific classes
            tf = false;
        end
        
        function result = diff_private(~,~)
            %default implementation for the derivative
            %the derivative of anything is a priori zero.
            %this function is overloaded in specific classes
            result = SYM_ZERO;
        end
        
        function Ta = get_factors_private(arg)
            Ta = arg;
        end
        function Ta = get_terms_private(arg)
            Ta = arg;
        end
        
        
        
        function result = inv_private(a)
            %default implementation for the inv of s scalar symbolic object
            result = simplify_private(SYM_FUNCTION('power',a,-SYM_ONE));
        end
        function result = simplify_private(a)
            %default implementation for simplification
            result = a;
        end
        function result = test_fct_private(a) %#ok<MANU>
            %default implementation for test_function
            %a variation of anything is a priori zero.
            %this function is overloaded is specific classes
            result = SYM_ZERO;
        end
        
        function [scalar, entity] = normalize_private(arg)
            %normalizes the argument to a scalar factor and a normalized
            %entity.
            %default implementation for normalize
            %this function can be overloaded in derived classes
            if arg.num
                scalar = arg.value;
                entity = SYM_ONE;
            else
                scalar = 1;
                entity = arg;
            end
        end
        
        function result = remove_products_of_unknown_fields_private(arg)
            result = arg;
        end
        
        function result = subs_private(arg,old,new)
            %default implementation for subs
            if arg == old
                result = new;
            else
                result = arg;
            end
        end
        function result = simplify_products_private(arg,~)
            result = arg;
        end
        function result = sym2double_private(arg)
            result = NaN(size(arg));
        end
        
    end
    
    methods (Static, Sealed, Access = protected)
        %This is to specify that unspecified elements in an heterogeneous
        %array of objects deriving from SYM_BASE should be initialized as
        %SYM_ZERO
        function default_object = getDefaultScalarElement
            default_object = SYM_ZERO;
        end
        
        function converted_object = convertObject(SYM_BASE,objectToConvert) %#ok<INUSL>
            if(isnumeric(objectToConvert))
                converted_object = conv2SYM(objectToConvert);
            else
                error('Cannot convert to SYM_BASE');
            end
        end
        
    end
    
    
end

function [a,b] = match_scalar(a,b)
%if one of the inputs is scalar the function uses repmat to make it
%match the size of the other input.
if isscalar(a)
    a = repmat(a,size(b));
    return;
end
if isscalar(b)
    b = repmat(b,size(a));
    return
end
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