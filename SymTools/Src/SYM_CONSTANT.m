%SYMBOLIC PDE TOOLBOX CLASS FOR REPRESENTING A SCALAR CONSTANT
%
%Class constructor
%Syntax:
% a = SYM_CONSTANT(S) : (S is a character array) symbolic
% representation of the named scalar constant S
%
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Adrien Leygue (Adrien.Leygue@ec-nantes.fr)
%
%SYM CLASSES PACKAGE
%
classdef SYM_CONSTANT < SYM_BASE
    %CLASS FOR THE SYMBOLIC REPRESENTATION OF SCALAR CONSTANTS
    properties
        id = '';
    end
    methods
        function obj = SYM_CONSTANT(id)
            %Class constructor
            %Syntax:
            % a = SYM_CONSTANT(S) : (S is a charracter array) symbolic
            % representation of the named constant S
            
            obj.name = 'Constant';
            if nargin==1
                obj.id = id;
            else
                error('SYM_CONSTANT accepts one and only one argument');
            end
        end
    end
    methods (Access = protected)
        %FOR COMMENTS ON THE FUNCTIONS PURPOSE, SEE THE SYM_BASE CLASS
        function result = display_tex_private(arg)
            result = arg.id;
        end
        function result = sym2str_private(arg)
            result = arg.id;
        end
        function result = eq_private(a,b)
            result = strcmp(a.id,b.id);
        end
        function result = export_private(arg)
            result = {arg.name arg.id};
        end
        function result = separate_private(arg,pattern)
            nd = numel(pattern);
            result = [arg ; SYM_ONE(nd,1)];
        end        
    end
end