%SYMBOLIC PDE TOOLBOX CLASS FOR REPRESENTING FUNCTIONS (SIN, COS,...)
%
%Class constructor
%Syntax:
% This class should be accessed through the overloaded operators
% obj = SYM_FUNCTION(name,argument,parameters)
% examples: SYM_FUNCTION('sin',x)
%           SYM_FUNCTION('power',x,y)
% currently implemented: sin, cos, exp, power
%
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Adrien Leygue (Adrien.Leygue@ec-nantes.fr)
%
%SYM CLASSES PACKAGE
%
classdef SYM_FUNCTION < SYM_BASE
    %CLASS FOR THE SYMBOLIC REPRESENTATION OF KNOWN SCALAR FUNCTIONS OF ONE OR SEVERAL
    %VARIABLES
    properties
        id = '';
        param = SYM_ZERO.empty;
        argument = SYM_ZERO.empty;
    end
    methods
        function obj = SYM_FUNCTION(id,argument,param)
            %Class constructor for the symbolic representation of
            %known functions
            %Syntax:
            % a = SYM_FUNCTION(name,argument): e.g. F = SYM_FIELD('cos',x)
            % name: a string identifying the function in the catalog
            % argument: a valid symbolic expression
            
            if nargin<2
                error('SYM_FUNCTION needs at least 2 inputs');
            end
            obj.name = 'Function';
            obj.id = id;
            obj.argument = argument;
            switch obj.id
                case 'sin'
                    obj.param = SYM_ZERO.empty;
                case 'cos'
                    obj.param = SYM_ZERO.empty;
                case 'exp'
                    obj.param = SYM_ZERO.empty;
                case 'power'
                    obj.param = param;
                    if isa(argument,'SYM_FUNCTION')
                        if strcmp(argument.id,'power')
                            obj.param = obj.param*argument.param;
                            if obj.param.num
                                obj.param = conv2SYM(sym2double(obj.param));
                            end
                            obj.argument = argument.argument;
                        end
                    end
                otherwise
                    error('unknown function');
            end
        end
    end
    methods (Access = protected)
        %FOR COMMENTS ON THE FUNCTIONS PURPOSE, SEE THE SYM_BASE CLASS
        
        function result = display_tex_private(arg)
            
            result = ['\mathrm{' arg.id '}'];
            %append argument
            result = [result '(' display_tex_private(arg.argument(1))];
            for i=2:numel(arg.argument)
                result = [result ',' display_tex_private(arg.argument(i))]; %#ok<AGROW>
            end
            for i=1:numel(arg.param)
                result = [result ',' display_tex_private(arg.param(i))]; %#ok<AGROW>
            end
            
            result = [result ')'];
        end
        
        
        function result = eq_private(a,b)
            %if the number of variables are different then false
            if (numel(a.argument) ~= numel(b.argument))
                result = false;
                return;
            end
            %if the number of parameters are different then false
            if (numel(a.argument) ~= numel(b.argument))
                result = false;
                return;
            end
            
            result =  all(a.argument==b.argument) & all(a.param==b.param);
        end
        
        function result = export_private(arg)
            result = {arg.name arg.id arg.argument.export arg.param.export};
        end
        
        %THESE FUNCTIONS ARE OVERLOADED FROM THE BASE CLASS
        
        function result = diff_private(arg,var)
            switch arg.id
                case 'sin'
                    result = SYM_FUNCTION('cos',arg.argument)*arg.argument.diff(var);
                case 'cos'
                    result = -SYM_FUNCTION('sin',arg.argument)*arg.argument.diff(var);
                case 'exp'
                    result = SYM_FUNCTION('exp',arg.argument)*arg.argument.diff(var);
                case 'power'
                    switch arg.param
                        case SYM_ZERO
                            result = SYM_ZERO;
                        case SYM_ONE
                            result = arg.argument.diff(var);
                        otherwise
                            result = arg.param*SYM_FUNCTION('power',arg.argument,arg.param-1)*arg.argument.diff(var);
                    end
                otherwise
                    error('unknown function');
            end
        end
        
        
        function result = separate_private(arg,pattern) %#ok<STOUT,INUSD>
            error('Cannot separate an arbitrary function');
        end
        
        function result = simplify_private(arg)
            result =arg;
            switch arg.id
                case 'sin'
                    if arg.argument==SYM_ZERO
                        result = SYM_ZERO;
                    end
                case 'cos'
                    if arg.argument==SYM_ZERO
                        result = SYM_ONE;
                    end
                case 'exp'
                    if arg.argument==SYM_ZERO
                        result = SYM_ONE;
                    end
                case 'power'
                    switch arg.param
                        case SYM_ZERO
                            result = SYM_ONE;
                        case SYM_ONE
                            result = arg.argument;
                    end
                    if (arg.argument.num && arg.param.num)
                        result = SYM_REAL(power(sym2double(arg.argument),sym2double(arg.param)));
                    end
                otherwise
                    error('unknown function');
            end
        end
        
        function result = subs_private(arg,old,new)
            %if the subs does not apply to the object itself, propagates to
            %the argument & parameter
            if arg == old
                result = new;
            else
                new_argument = subs(arg.argument,old,new);
                if ~isempty(arg.param)
                    new_param = subs(arg.param,old,new);
                    result = simplify_private(SYM_FUNCTION(arg.id,new_argument,new_param));
                else
                    result = simplify_private(SYM_FUNCTION(arg.id,new_argument));
                end
            end
            
        end
        
        function result = simplify_products_private(arg,opt)
            if isempty(arg.param)
                result = SYM_FUNCTION(arg.id,simplify_products(arg.argument,opt));
            else
                result = SYM_FUNCTION(arg.id,simplify_products(arg.argument,opt),simplify_products(arg.param,opt));
            end
        end
        
        function result = sym2str_private(arg)
            result = arg.id;
            %append argument
            result = [result '(' sym2str(arg.argument(1))];
            for i=2:numel(arg.argument)
                result = [result ',' sym2str(arg.argument(i))]; %#ok<AGROW>
            end
            for i=1:numel(arg.param)
                result = [result ',' sym2str(arg.param(i))]; %#ok<AGROW>
            end
            
            result = [result ')'];
        end

    end
end