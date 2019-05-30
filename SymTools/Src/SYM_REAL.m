%SYMBOLIC PDE TOOLBOX CLASS FOR REPRESENTING REALS
%
%Class constructor
%Syntax:
% a = SYM_REAL  : equivalent to SYM_REAL(0)
% a = SYM_REAL(V) : (V: real scalar) symbolic representation of V
% a = SYM_REAL(V) : (V: real array) array of SYM_REAL of the
% same size as V and such that a(i,j,..,k) = SYM_REAL(V(i,j,...,k))
%
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Adrien Leygue (Adrien.Leygue@ec-nantes.fr)
%
%SYM CLASSES PACKAGE
%
classdef SYM_REAL < SYM_BASE
    %CLASS FOR THE SYMBOLIC REPRESENTATION OF ZEROS
    properties
        value = 0;
    end
    methods
        function obj = SYM_REAL(arg)
            %Class constructor
            %Syntax:
            % a = SYM_REAL  : equivalent to SYM_REAL(0)
            % a = SYM_REAL(V) : (V: real scalar) symbolic representation of V
            % a = SYM_REAL(V) : (V: real array) array of SYM_REAL of the
            % same size as V and such that a(i,j,..,k) = SYM_REAL(V(i,j,...,k))
            
            %no argument
            if nargin==0
                obj.name = 'Real';
                obj.num = true;
                return;
            end
            %one argument
            if isscalar(arg)
                obj.name = 'Real';
                obj.num = true;
                obj.value = arg;
                return;
            else
                obj(numel(arg)) = SYM_REAL(arg(end));
                for i=1:numel(arg)
                    obj(i) = SYM_REAL(arg(i));
                end
                obj = reshape(obj,size(arg));
                return;
            end
        end
    end
    methods (Access = protected)
        %FOR COMMENTS ON THE FUNCTIONS PURPOSE, SEE THE SYM_BASE CLASS
        function result = display_tex_private(arg)
            result = num2str(arg.value,'%e');
        end
        
        function result = eq_private(a,b)
            result = (a.value==b.value);
        end
        function result = export_private(arg)
            result = {arg.name arg.value};
        end
        function result = inv_private(arg)
            result = SYM_REAL(arg.value.^(-1));
        end
        function result = separate_private(arg,pattern)
            nd = numel(pattern);
            result = [arg; SYM_ONE(nd,1)];
        end
        function result = sym2double_private(arg)
            result = [arg.value];
        end
        function result = sym2str_private(arg)
            result = num2str(arg.value,'%e');
        end
        
        
    end
end