%SYMBOLIC PDE TOOLBOX CLASS FOR REPRESENTING a SUM
%Syntax:
% a = SYM_SUM(b,c): creates the symbolic representation of b+c
% NB: inputs cannot be of type SYM_ZERO (this is checked in atomic_plus)
% NB: call collect_scalars to simplify the sum
%
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Adrien Leygue (Adrien.Leygue@ec-nantes.fr)
%

classdef SYM_SUM < SYM_BASE
    properties
        Nterms = 0;
        terms = SYM_ZERO.empty;
    end
    methods
        function obj = SYM_SUM(a,b)
            %Class constructor for the symbolic representation of sums
            %This constructor is only called through the overloaded plus,
            %mtimes,..operator or through the atomic_plus function
            %Syntax:
            % a = SYM_SUM(b,c): creates the symbolic representation of b+c
            % NB: inputs cannot be of type SYM_ZERO (this is checked in atomic_plus)
            % NB: call collect_scalars to simplify the sum
            obj.name = 'Sum';
            if nargin==0
                return;
            end
            
            assert(nargin==2,'SYM_SUM constructor accepts 0 or 2 arguments');
            
            %in case a or b is a sum we extract the terms
            Ta = get_terms_private(a);
            Tb = get_terms_private(b);
            %concatenate all
            obj.terms = [Ta Tb];
            obj.Nterms = numel(obj.terms);
        end

    end
    methods (Access = protected)
        %FOR COMMENTS ON THE FUNCTIONS PURPOSE, SEE THE SYM_BASE CLASS
        function result = display_tex_private(arg)
            result = display_tex_private(arg.terms(1));
            for i=2:arg.Nterms
                result = [result ' + ' display_tex_private(arg.terms(i))]; %#ok<AGROW>
            end
        end
        
        function result = eq_private(a,b)
            if (a.Nterms ~=b.Nterms)
                result = false;
                return;
            end
            result = all(a.terms == b.terms);
        end
        function result = export_private(arg)
            result = {arg.name export(arg.terms)};
        end
        function result = separate_private(arg,pattern)
            result = SYM_ZERO(numel(pattern)+1,arg.Nterms);
            for i=1:arg.Nterms
                result(:,i) = separate_private(arg.terms(i),pattern);
            end
%            error('One should only separate individual terms');
        end
        
        %THESE FUNCTIONS ARE OVERLOADED FROM THE BASE CLASS
        function tf = contains_tagged_field_private(arg,tag)
            %propagates the query to the terms
            tmp = false(size(arg.terms));
            for i=1:arg.Nterms
                tmp(i) = contains_tagged_field_private(arg.terms(i),tag);
            end
            tf = any(tmp);
        end
        function result = diff_private(arg,var)
            %diff of sum is equal to the sum of diff
            result = SYM_ZERO;
            for i=1:arg.Nterms
                result = atomic_plus(result,diff_private(arg.terms(i),var));
            end
            if isa(result,'SYM_SUM')
                result = collect_scalars(result);
            end
        end
        function result = test_fct_private(arg)
            %test_fct od a sum is the sum of the test_fct (similar to diff)
            result = SYM_ZERO;
            for i=1:arg.Nterms
                result = atomic_plus(result,test_fct_private(arg.terms(i)));
            end
            if isa(result,'SYM_SUM')
                result = collect_scalars(result);
            end
        end
        function result = subs_private(arg,old,new)
            %if the subs does not apply to the object itself, propagates to
            %the terms
            if arg == old
                result = new;
            else
                result = SYM_ZERO;
                for i=1:arg.Nterms
                    result = atomic_plus(result,subs_private(arg.terms(i),old,new));
                end
                if isa(result,'SYM_SUM')
                    result = collect_scalars(result);
                end
            end
        end
        
        function result = simplify_products_private(arg,opt)
                result = SYM_ZERO;
                for i=1:arg.Nterms
                    result = atomic_plus(result,simplify_products_private(arg.terms(i),opt));
                end
                if isa(result,'SYM_SUM')
                    result = collect_scalars(result);
                end
        end
        
        function result = remove_products_of_unknown_fields_private(arg)
                result = SYM_ZERO;
                for i=1:arg.Nterms
                    result = atomic_plus(result,remove_products_of_unknown_fields_private(arg.terms(i)));
                end
                if isa(result,'SYM_SUM')
                    result = collect_scalars(result);
                end
        end
        
        function result = sym2str_private(arg)
            result = sym2str(arg.terms(1));
            for i=2:arg.Nterms
                result = [result '+' sym2str(arg.terms(i))]; %#ok<AGROW>
            end
        end
        function Ta = get_terms_private(arg)
            Ta = arg.terms;
        end
    end
    methods
        function result = collect_scalars(arg)
            %Identifies terms that differ only by a scalar factor and
            %regtoups them
            %normalization
            [tmp_scalars,tmp_entities] = normalize(arg.terms);
            %ordering
            [sorted_entities,str,p] = sort_disp(tmp_entities);
            sorted_scalars = tmp_scalars(p);
            
            %identification of groups
            [~,gather_vect,scatter_vect] = unique(str,'stable');
            %new terms are computed
            new_terms = sorted_entities(gather_vect);
            new_scalars = accumarray(scatter_vect,sorted_scalars,size(gather_vect));
            new_terms = times(new_terms,new_scalars');
            
            mask =true(size(new_terms));
            for i=1:numel(new_terms)
                mask(i) = ~isa(new_terms(i),'SYM_ZERO');
            end
            new_terms = new_terms(mask);
            %the output is constructed
            result = SYM_SUM;
            result.terms = new_terms;
            result.Nterms = numel(new_terms);
            switch result.Nterms
                case 1
                result = result.terms;
                case 0
                    result = SYM_ZERO;
            end
        end
    end
    
end



