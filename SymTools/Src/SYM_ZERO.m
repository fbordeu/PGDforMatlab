%SYMBOLIC PDE TOOLBOX CLASS FOR REPRESENTING ZERO
%
%Class constructor
%Syntax:
% a = SYM_ZERO  : a single object representing ZERO
% a = SYM_ZERO(M) : (scalar M )an M by M array of SYM_ZERO
% a = SYM_ZERO(S) : (vector S) an array of SYM_ZERO with S = size(a)
% a = SYM_ZERO(M,N) : an M by N array of SYM_ZERO
% a = SYM_ZERO(M,N,...P) : an M by N by... by P array of SYM_ZERO
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
classdef SYM_ZERO < SYM_BASE
    %CLASS FOR THE SYMBOLIC REPRESENTATION OF ZEROS
    properties
        value = 0;
    end
    methods
        function obj = SYM_ZERO(varargin)
            %Class constructor
            %Syntax:
            % a = SYM_ZERO  : a single object
            % a = SYM_ZERO(M) : (scalar M )an M by M array of SYM_ZERO
            % a = SYM_ZERO(S) : (vector S) an array of SYM_ZERO with S = size(a)
            % a = SYM_ZERO(M,N) : an M by N array of SYM_ZERO
            % a = SYM_ZERO(M,N,...P) : an M by N by... by P array of SYM_ZERO
            
            %no argument
            if isempty(varargin)
                %init the inherited properties
                obj.name = 'Zero';
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
                    error('SYM_ZERO with one input requires a scalar or vector')
                end
                %more than one argument
            else
                sz = [varargin{:}];
            end
            
            obj(prod(sz)) = SYM_ZERO; %OK as SYM_ZERO is the default scalar element
            obj = reshape(obj,sz);
        end
        
    end
    methods (Access = protected)
        %FOR COMMENTS ON THE FUNCTIONS PURPOSE, SEE THE SYM_BASE CLASS
        function result = display_tex_private(~)
            result = '0';
        end
        function result = eq_private(~,~)
            result = true;
        end
        function result = export_private(arg)
            result = {arg.name};
        end
        function result = inv_private(~) %#ok<STOUT>
            error('SYM_ZERO: Error cannot compute inv(SYM_ZERO)');
        end
        function result = separate_private(~,pattern)
            nd = numel(pattern);
            result = SYM_ZERO(nd+1,1);
        end
        function result = sym2double_private(arg)
            result = [arg.value];
        end
        function result = sym2str_private(~)
            result = '0';
        end
    end
end