%SYMBOLIC PDE TOOLBOX CLASS FOR REPRESENTING A PRODUCT
%
%Class constructor
%Syntax:
% a = SYM_PROD(b,c): creates the symbolic representation of b*c
% NB: inputs cannot be of type SYM_ZERO or SYM_ONE
%
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Adrien Leygue (Adrien.Leygue@ec-nantes.fr)
%
%SYM CLASSES PACKAGE
%
%
%Last modification of this file: 06 sept. 2013
classdef SYM_PROD < SYM_BASE
    properties
        Nfactors = 0;
        factors = SYM_ZERO.empty;
    end
    methods
        function obj = SYM_PROD(a,b)
            %Class constructor for the symbolic representation of a product
            %This constructor is only called through the overloaded times,
            %mtimes,... operators or through the atomic_plus function
            %Syntax:
            % a = SYM_PROD(b,c): creates the symbolic representation of b*c
            % NB: inputs cannot be of type SYM_ZERO or SYM_ONE (this is checked in atomic_plus)
            obj.name = 'Prod';
            
            if nargin==0
                return;
            end
            assert(nargin==2,'SYM_PROD constructor accepts 0 or 2 arguments');
            
            Fa = get_factors_private(a);
            Fb = get_factors_private(b);
            
            obj.factors = [Fa Fb];
            obj.Nfactors = numel(obj.factors);
            
            mask = [obj.factors.num];
            if nnz(mask)>1
                obj.factors =[SYM_REAL(prod([obj.factors(mask).value])) obj.factors(~mask)];
            end
            obj.factors = sort_disp(obj.factors);
            obj.Nfactors = numel(obj.factors);
        end
    end
    methods (Access = protected)
        %FOR COMMENTS ON THE FUNCTIONS PURPOSE, SEE THE SYM_BASE CLASS
        function result = display_tex_private(arg)
            result = display_tex_private(arg.factors(1));
            for i=2:arg.Nfactors
                result = [result ' \cdot ' display_tex_private(arg.factors(i))]; %#ok<AGROW>
            end
        end
        function result = eq_private(a,b)
            if (a.Nfactors ~=b.Nfactors)
                result = false;
                return;
            end
            result = all(a.factors == b.factors);
        end
        function result = export_private(arg)
            result = {arg.name export(arg.factors)};
        end
        %THESE FUNCTIONS ARE OVERLOADED FROM THE BASE CLASS
        function [scalar, entity] = normalize_private(arg)
            %normalizes each factor then multiplies everything together
            [tmp_scalars, tmp_entities] = normalize(arg.factors);
            scalar = prod(tmp_scalars);
            entity = SYM_ONE;
            for i=1:arg.Nfactors
                entity = atomic_times(entity,tmp_entities(i));
            end
        end
        
        function result = diff_private(arg,var)
            %applies the rule for the derivative of a product
            [a,b] = split_prod(arg);
            result = plus(atomic_times(diff_private(a,var),b),atomic_times(diff_private(b,var),a)) ;
        end
        function result = test_fct_private(arg)
            %test_function of a product (identical to derivative)
            [a,b] = split_prod(arg);
            result = plus(atomic_times(test_fct_private(a),b),atomic_times(test_fct_private(b),a)) ;
        end
        
        function result = subs_private(arg,old,new)
            %if the subs does not apply to the object itself, propagates to
            %the factors
            if arg == old
                result = new;
            else
                result = SYM_ONE;
                for i=1:arg.Nfactors
                    result = atomic_times(result,subs_private(arg.factors(i),old,new));
                end
            end
        end
        
        function result = remove_products_of_unknown_fields_private(arg)
            
            tmp = false(size(arg.factors));
            for i=1:arg.Nfactors
                tmp(i) = contains_tagged_field_private(arg.factors(i),'unknown');
            end
            
            if nnz(tmp)>1
                 result=SYM_ZERO();
             else
                 result = arg;
             end
        end
        
        function result = remove_SYM_ONES(arg)
            mask = ~(arg.factors==SYM_ONE);
            if nnz(mask)==1
                result=arg(mask);
            else
                result = arg;
                result.factors = result.factors(mask);
                result.Nfactors = nnz(mask);
            end
        end
        
        function result = simplify_products_private(arg,opt)
            tmp_exponents = SYM_ONE(arg.Nfactors);
            tmp_entities = arg.factors;
            for f = 1:arg.Nfactors
                if isa(tmp_entities(f),'SYM_FUNCTION') && strcmp(tmp_entities(f).id,'power')
                    tmp_exponents(f) = tmp_entities(f).param;
                    tmp_entities(f) = tmp_entities(f).argument;
                end
            end
            [sorted_entities,str,p] = sort_disp(tmp_entities);
            sorted_exponents = tmp_exponents(p);
            %identification of groups
            [~,gather_vect,scatter_vect] = unique(str,'stable');
            new_factors = sorted_entities(gather_vect);
            new_exponents = SYM_ZERO(size(new_factors));
            for i=1:numel(sorted_exponents)
                new_exponents(scatter_vect(i)) = new_exponents(scatter_vect(i)) + sorted_exponents(i);
            end

            mask =true(size(new_exponents));
            for f=1:numel(new_exponents)
                mask(f) = ~isa(new_exponents(f),'SYM_ZERO');
            end
            result = SYM_PROD;
            
            for f = find(mask)
                    switch opt
                        case 'collect'
                            if isa(new_exponents(f),'SYM_ONE')
                                result.factors(result.Nfactors+1) = new_factors(f);
                            else 
                            result.factors(result.Nfactors+1) = SYM_FUNCTION('power',new_factors(f),new_exponents(f));
                            end
                            result.Nfactors = result.Nfactors+1;
                        case 'expand'
                            if (new_exponents(f).num) && (new_exponents(f).value==fix(new_exponents(f).value)) && (new_exponents(f).value>0)
                                result.factors(result.Nfactors+(1:new_exponents(f).value)) = new_factors(f);
                                result.Nfactors = result.Nfactors+new_exponents(f).value;
                            else 
                            result.factors(result.Nfactors+1) = SYM_FUNCTION('power',new_factors(f),new_exponents(f));
                            result.Nfactors = result.Nfactors+1;
                            end
                        otherwise
                            error('unknown option for simplify_products');
                    end
            end
            if result.Nfactors==0
                result = SYM_ONE;
            end
        end
        
        function result = separate_private(arg,pattern)
            %computes the product (times) of the separation of each factor
            nd = numel(pattern);
            result = SYM_ONE(nd+1,1);
            for i=1:arg.Nfactors
                result = result.*separate_private(arg.factors(i),pattern);
            end
        end
        function tf = contains_tagged_field_private(arg,tag)
            %propagates the query to the factors
            tmp = false(size(arg.factors));
            for i=1:arg.Nfactors
                tmp(i) = contains_tagged_field_private(arg.factors(i),tag);
            end
            tf = any(tmp);
        end
        function result = sym2str_private(arg)
            result = sym2str_private(arg.factors(1));
            for i=2:arg.Nfactors
                result = [result '*' sym2str_private(arg.factors(i))]; %#ok<AGROW>
            end
        end
        function Ta = get_factors_private(arg)
                Ta = arg.factors;
        end


    end
end


function [a,b] = split_prod(arg)
    %decomposes a product in the first factor (a) times the product of the
    %other factors (b) 
    a = arg.factors(1);
    if arg.Nfactors == 2
        b = arg.factors(2);
    else
        b = arg;
        b.factors = b.factors(2:end);
        b.Nfactors = b.Nfactors-1;
    end
end


