%SYMBOLIC PDE TOOLBOX CLASS FOR REPRESENTING A FIELD
%
%Class constructor
%Syntax:
% a = SYM_FIELD(name,variables,tag,Ncomp,comp): e.g. F = SYM_FIELD('F',[x y],'unknown',3,2)
% name: a string identifying the field
% variables: an array of SYM_COORD of the variables of the field
% tag: a string specifying the nature of the field:
%      - ''              : a known function
%      - 'unknown' or 'u': an unknown function
%      - 'test'    or 't': a test function in a weak formulation
% Ncomp: number of components of the field (e.g. 3 for a velocity field)
% comp : specific component of the field to build (comp out of Ncomp). If
% ommitted, comp=1:Ncomp
%
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Adrien Leygue (Adrien.Leygue@ec-nantes.fr)
%
%
classdef SYM_FIELD < SYM_BASE
    %CLASS FOR THE SYMBOLIC REPRESENTATION OF SCALAR FUNCTIONS/FIELDS OF ONE OR SEVERAL
    %VARIABLES
    properties
        variables = SYM_ZERO.empty;
        id = '';
        diff_symbol = [];
        depend_symbol = [];
        tag = '';
        comp = 1;
        Ncomp = 1;
        str = '';
    end
    methods
        function obj = SYM_FIELD(id,variables,tag,Ncomp,comp)
            
            if nargin==0
                error('SYM_FIELD needs at least 2 inputs');
            end
            obj.name = 'Field';
            
            if nargin<3
                tag = 'known';
            end
            
            if nargin<4
                Ncomp = 1;
            end
            
            obj.id = id;
            obj.variables = variables;
            obj.depend_symbol = ones(size(variables));
            obj.diff_symbol = zeros(1,numel(obj.variables));
            obj.Ncomp = Ncomp;
            
            if numel(tag)==1
                switch tag
                    case 'u'
                        tag = 'unknown';
                    case 't'
                        tag = 'test';
                    case 'k'
                        tag = 'known';
                end
            end
            if isempty(tag)
                tag = 'known';
            end
            
            assert(ismember(tag,{'','known','test','unknown'}),'SYM_FIELD: thirt input is a string among: '''' , ''known'', ''test'' or ''unknown'' ' );
            obj.tag = tag;
            obj = repmat(obj,[Ncomp 1]);
            for i=1:Ncomp
                obj(i).comp = i;
            end
            if(nargin>4)
                assert(max(comp)<=Ncomp);
                obj = obj(comp);
            end
            for i=1:numel(obj)
                obj(i).str = create_str(obj(i));
            end
        end
        function result = set_dependency(arg,dvar)
            result = arg;
            for k=1:numel(arg)
            result(k).depend_symbol = zeros(size(result(k).depend_symbol));
            for i=1:numel(dvar)
                mask = (dvar(i)==result(k).variables);
                if ~any(mask)
                    error('incorrect input variable');
                else
                    result(k).depend_symbol(mask) = 1;
                end
            end
            end
        end
    end
    methods (Access = protected)
        %FOR COMMENTS ON THE FUNCTIONS PURPOSE, SEE THE SYM_BASE CLASS
        
        function result = display_tex_private(arg)
            
            diff_str = '';
            %generate the string for the derivative
            if any(arg.diff_symbol)
                %numerator
                if sum(arg.diff_symbol)>1
                    diff_str = ['\frac{\partial^' num2str(sum(arg.diff_symbol)) '}{'];
                else
                    diff_str = '\frac{\partial}{';
                end
                %denominator
                for i=1:numel(arg.variables)
                    if arg.diff_symbol(i)>1
                        diff_str = [diff_str '\partial ' arg.variables(i).display_tex_private '^' num2str(arg.diff_symbol(i))]; %#ok<AGROW>
                    elseif arg.diff_symbol(i)==1
                        diff_str = [diff_str ' \partial ' arg.variables(i).display_tex_private]; %#ok<AGROW>
                    end
                end
                diff_str = [diff_str '} '];
            end
            result = [diff_str arg.id];
            if arg.Ncomp>1
                result = [result '[' num2str(arg.comp) ']'];
            end
            %modify if test function
            if strcmpi(arg.tag,'test')
                result = [result '^*'];
            end
            %append variables
            result = [result '(' arg.variables(1).display_tex_private];
            for i=2:numel(arg.variables)
                result = [result ',' arg.variables(i).display_tex_private]; %#ok<AGROW>
            end
            result = [result ')'];
        end
        
        
        
        function result = eq_private(a,b)
            %if the number of variables are different then false
            if (numel(a.variables) ~= numel(b.variables))
                result = false;
                return;
            end
            %variables can be compared as the cell array are of the same
            %size
            result = strcmp(a.id,b.id)&strcmp(a.tag,b.tag)&(a.comp==b.comp)&(a.Ncomp==b.Ncomp)&all(a.diff_symbol==b.diff_symbol)& all(a.variables==b.variables);
        end
        
        function result = export_private(arg)
            result = {arg.name arg.id arg.variables.export arg.diff_symbol arg.tag arg.comp arg.Ncomp};
        end
        
        %THESE FUNCTIONS ARE OVERLOADED FROM THE BASE CLASS
        
        function result = diff_private(arg,var)
            %update the diff_symbol
            idx = (var==arg.variables);
            if any(idx) && arg.depend_symbol(idx)
                result = arg;
                result.diff_symbol(idx) = result.diff_symbol(idx)+1;
                result.str = create_str(result);
            else
                result = SYM_ZERO;
            end
        end
        function result = test_fct_private(arg)
            %a known function becomes SYM_ZERO
            if strcmpi(arg.tag,'known')
                result = SYM_ZERO;
                return;
            end
            %an unknown function becomes a test function
            if strcmpi(arg.tag,'unknown')
                result = arg;
                result.tag = 'test';
                result.str = create_str(result);
                return;
            end
            %error if we already have a test_function
            error(['test_fct is inappropriate on a field tagged ' arg.tag]);
        end
        
        function tf = contains_tagged_field_private(arg,tag)
            %test the tag
            tf = strcmp(arg.tag,tag);
        end
        
        function result = separate_private(arg,pattern)
            nd = numel(pattern);
            sep_vars = [];
            for i=1:nd
                sep_vars = [sep_vars pattern{i}]; %#ok<AGROW>
            end
            sep_vars = sort_disp(sep_vars);
            assert(all(sep_vars==sort_disp(arg.variables)),'Variables mismatch in the separation');
            
            result = SYM_ONE(nd+1,1);
            for i=1:nd
                %for each group of variables we identify they position in arg.variables
                id_pattern = zeros(size(pattern{i}));
                for j=1:numel(pattern{i})
                    pos = find(pattern{i}(j)==arg.variables);
                    id_pattern(j) = pos;
                end
                id_pattern = sort(id_pattern);
                result(i+1) = arg;
                result(i+1).variables = arg.variables(id_pattern);
                result(i+1).diff_symbol = arg.diff_symbol(id_pattern);
                result(i+1).depend_symbol = arg.depend_symbol(id_pattern);
                result(i+1).str = create_str(result(i+1));
            end
        end
        
        function result = sym2str_private(arg)
            result = arg.str;
        end
        
        function result = create_str(arg)
            %function base name
            tmp = '';
            if arg.Ncomp>1
                tmp = ['[' num2str(arg.comp) ']'];
            end
            result = [arg.id tmp '(' arg.variables(1).sym2str];
            %append variables
            for i=2:numel(arg.variables)
                result = [result ',' arg.variables(i).sym2str]; %#ok<AGROW>
            end
            result = [result ')'];
            %modify if test function
            if strcmpi(arg.tag,'test')
                result = [result '°'];
            end
            %append derivatives
            tmp = arg.diff_symbol;
            while(any(tmp))
                pos = find(tmp,1,'first');
                tmp(pos) = tmp(pos)-1;
                result = [result '_' arg.variables(pos).sym2str]; %#ok<AGROW>
            end
        end
    end
end