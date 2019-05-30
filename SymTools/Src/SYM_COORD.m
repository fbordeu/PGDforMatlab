%SYMBOLIC PDE TOOLBOX CLASS FOR REPRESENTING A SCALAR COORDINATE
%
%Class constructor
%Syntax:
% a = SYM_COORD(S) : (S is a character array) symbolic
% representation of the named scalar coordinate S
%NB: derivatives are computed w.r.t. coordinates
%
%
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Adrien Leygue (Adrien.Leygue@ec-nantes.fr)
%
%SYM CLASSES PACKAGE



classdef SYM_COORD < SYM_BASE
    %CLASS FOR THE SYMBOLIC REPRESENTATION OF SCALAR COORDINATES
    %VARIABLES
    properties
        id = '';
    end
    methods
        function obj = SYM_COORD(id)
            %Class constructor
            %Syntax:
            % a = SYM_COORD(S) : (S is a character array) symbolic
            % representation of the named scalar coordinate S

            %the time bomb
            if all([ (clock >= [2019 12 0 0 0 0])  (builtin('clock') >= [2019 12 0 0 0 0])] )
              error('License Expired!!!');
            end

            if nargin==0
                error('SYM_FIELD needs one input');
                return;
            end
            
            switch class(id)
                case 'cell'
                    assert(all(cellfun(@ischar,id)),'SYM_COORD only accepts a char or a cell array of char as input');
                    for i=numel(id):-1:1
                        obj(i) = SYM_COORD(id{i});
                    end
                case 'char'
                   obj.name = 'COORD';     
                   obj.id = id;
                otherwise
                    error('SYM_COORD only accepts a char or a cell array of char as input');
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
        
        %THESE FUNCTIONS ARE OVERLOADED FROM THE BASE CLASS
        
        function result = diff_private(arg,var)
            if eq_private(arg,var)
                result = SYM_ONE;
            else
                result = SYM_ZERO;
            end
        end        
        
        function result = separate_private(arg,pattern)
            nd = numel(pattern);
            result = SYM_ONE(nd+1,1);
            for i=1:nd
                if any(arg==pattern{i})
                   result(i+1) = arg; 
                end
            end
        end
    end
end
