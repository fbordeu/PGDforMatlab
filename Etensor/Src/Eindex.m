%Index class Tensor expression with Einstein notation for the indices
%see examples for usage.
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Adrien Leygue (Adrien.Leygue@ec-nantes.fr)
%

classdef Eindex < handle
    properties
        symbol
        range
        value
    end
    methods
        function obj = Eindex(symbol,range)
            narginchk(2,2);
            validateattributes(symbol,{'char'},{},mfilename,'symbol',1);
            validateattributes(range,{'numeric'},{'scalar','integer','finite','nonnegative'},mfilename,'range',2);
            obj.symbol = symbol;
            obj.range = range;
            obj.value = 1;
        end
        function disp(obj)
            a = '';
            for i=1:numel(obj)
                a = sprintf('%s %s[1:%d]',a,obj(i).symbol,obj(i).range);
            end
            disp(a);
        end
        function reset(obj)
            [obj.value] = deal(1);
        end
        function result = getRelevant(obj)
            result = (obj([obj.range]>1));
        end
        function result = ismax(obj)
            result = all([obj.value]==[obj.range]);
        end
        function increment(obj)
            assert(~ismax(obj),'cannot increment')
            if numel(obj)==1
                obj.value = obj.value+1;
                return;
            end
            sz = [obj.range];
            tmp = cell(1,numel(obj));
            [tmp{:}] = ind2sub(sz,sub2ind(sz,obj.value)+1);
            [obj.value] = deal(tmp{:});
        end
        function [result,factor] = voigtIndex(obj,sample1)
            %obj should be of length 2 or 4
            %range should be 3
            switch numel(obj)
                case 4
                    [result1,factor1] = voigtIndex(obj(1:2),sample1);
                    [result2,factor2] = voigtIndex(obj(3:4),sample1);
                    result = [result1(1) result2(1)];
                    factor = factor1*factor2;
                case 2
                    tmp = [obj.value];
                    if tmp(1)==tmp(2)
                        result = [tmp(1) 1];
                        factor = sample1;
                    elseif sum(tmp)==5
                        result= [4 1];
                        factor = sqrt(2*sample1);
                    elseif sum(tmp)==4
                        result= [5 1];
                        factor = sqrt(2*sample1);
                    elseif sum(tmp)==3
                        result = [6 1];
                        factor = sqrt(2*sample1);
                    else
                        error('surprising error');
                    end
                otherwise
                    error('voigtIndex only accepts a multiindex of length 2 or 4');
            end  
        end
        
    end
end