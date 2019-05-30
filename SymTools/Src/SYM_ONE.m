%SYMBOLIC PDE TOOLBOX CLASS FOR REPRESENTING ONE
%
%Class constructor
%Syntax:
% a = SYM_ONE  : a single object
% a = SYM_ONE(M) : (scalar M )an M by M array of SYM_ONE
% a = SYM_ONE(S) : (vector S) an array of SYM_ONE with S = size(a)
% a = SYM_ONE(M,N) : an M by N array of SYM_ONE
% a = SYM_ONE(M,N,...P) : an M by N by... by P array of SYM_ONE
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
%Last modification of this file: 24 oct. 2013
classdef SYM_ONE < SYM_BASE
    %CLASS FOR THE SYMBOLIC REPRESENTATION OF ONES
    properties
        value = 1;
    end
    methods
        function obj = SYM_ONE(varargin)
            %Class constructor
            %Syntax:
            % a = SYM_ONE  : a single object
            % a = SYM_ONE(M) : (scalar M )an M by M array of SYM_ONE
            % a = SYM_ONE(S) : (vector S) an array of SYM_ONE with S = size(a)
            % a = SYM_ONE(M,N) : an M by N array of SYM_ONE
            % a = SYM_ONE(M,N,...P) : an M by N by... by P array of SYM_ONE
            
            %no argument
            if isempty(varargin)
                %init the inherited properties
                obj.name = 'One';
                obj.num = true;
                return;
            end
            
            %one argument
            if numel(varargin)==1
                tmp = varargin{1};
                if isscalar(tmp)
                    sz = [tmp tmp];
                elseif isvector(varargin{1})
                    sz = tmp(:)';
                else
                    error('SYM_ONE with one input requires a scalar or vector')
                end
                %more than one argument
            else
                sz = [varargin{:}];
            end
            
            obj(prod(sz)) = SYM_ONE;
            %explicit initialization is needed as the default scalar
            %element is SYM_ZERO
            for i=1:numel(obj)
                obj(i) = SYM_ONE;
            end
            obj = reshape(obj,sz);
        end
        
    end
    methods (Access = protected)
        %FOR COMMENTS ON THE FUNCTIONS PURPOSE, SEE THE SYM_BASE CLASS
        function result = display_tex_private(~)
            result = '1';
        end
        
        function result = eq_private(~,~)
            result = true;
        end
        function result = export_private(arg)
            result = {arg.name};
        end
        function result = inv_private(arg)
            result = arg;
        end
        function result = separate_private(~,pattern)
            nd = numel(pattern);
            result = SYM_ONE(nd+1,1);
        end
        function result = sym2double_private(arg)
            result = [arg.value];
        end
        function result = sym2str_private(~)
            result = '1';
        end
        
    end
end