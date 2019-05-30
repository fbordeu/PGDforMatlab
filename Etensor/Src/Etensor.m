%Tensor expression with Einstein notation for the indices
%see examples for usage.
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Adrien Leygue (Adrien.Leygue@ec-nantes.fr)
%
classdef Etensor
    properties
        data
        sz
        indices
    end
    methods
        function obj = Etensor(data,varargin)
            narginchk(1,Inf);
            id = Eindex.empty;
            for i=1:numel(varargin)
                assert(isa(varargin{i},'Eindex'),'input 2 to end should be Eindex');
                id = [id;varargin{i}(:)]; %#ok<AGROW>
            end
            mask = [id.range]>1;
            switch nnz(mask)
                case 0
                    id = [Eindex('',1);Eindex('',1)];
                case 1
                    id = [id(mask);Eindex('',1)];
                otherwise
                    id = id(mask);
            end
            
            obj.indices = id(:);
            obj.sz = [id.range];
            assert(isequal(size(data),obj.sz),'inconsistent data size');
            obj.data = data;
        end
        function disp(obj)
            if numel(obj)>1
                for i=1:numel(obj)
                    obj(i).disp();
                end
                return;
            end
            id = obj.indices(obj.sz>1);
            a = sprintf('Etensor[%s]\n',strjoin({id.symbol},','));            
            disp(a);
            disp(obj.data);
        end
        function result = idx(obj,varargin)
            if numel(obj)>1
                for i=1:numel(obj)
                    result(i) = obj(i).idx(varargin{:}); %#ok<AGROW>
                end
                return;
            end
            id = Eindex.empty;
            for i=1:numel(varargin)
                assert(isa(varargin{i},'Eindex'),'input 2 to end should be Eindex');
                id = [id;varargin{i}(:)]; %#ok<AGROW>
            end
            id = id.getRelevant();
            old_id = obj.indices.getRelevant();
            
            assert(isequal([id.range],[old_id.range]),'Inconsistent indices range');
            result = Etensor(obj.data,id);
        end
        function result = permute(obj,varargin)
            if numel(obj)>1
                for i=1:numel(obj)
                    result(i) = obj(i).permute(varargin{:}); %#ok<AGROW>
                end
                return;
            end
            id = Eindex.empty;
            for i=1:numel(varargin)
                assert(isa(varargin{i},'Eindex'),'input 2 to end should be Eindex');
                id = [id;varargin{i}(:)]; %#ok<AGROW>
            end
            switch nnz(obj.sz>1)
                case 0
                    assert(isempty(id.getRelevant()),'inconsistent indices');
                    result = obj;
                case 1
                    assert(isequal(id.getRelevant(),obj.indices.getRelevant()),'inconsistent indices');
                    result = obj;
                otherwise
                    assert(isequal(sort(id),sort(obj.indices)),'inconsistent indices');
                    [~,locb] = ismember(id,obj.indices);
                    result = obj;
                    result.indices = id;
                    result.sz = obj.sz(locb);
                    result.data = permute(obj.data,locb);
            end
        end
        
        function result = contraction(varargin)
            T = Etensor.empty;
            for i=1:numel(varargin)
                assert(isa(varargin{i},'Etensor'),'All inputs should be Etensors');
                T = [T;varargin{i}]; %#ok<AGROW>
            end
            [all_idx,gather_vect,scatter_vect] = unique(vertcat(T.indices),'stable');
            occurences = accumarray(scatter_vect,1,size(gather_vect));
            mute_indices = all_idx(occurences>1);
            rem_indices = all_idx(occurences==1);
            mask = ([rem_indices.range]==1);
            mute_indices = [mute_indices;rem_indices(mask)];
            rem_indices = rem_indices(~mask);
            switch(numel(rem_indices))
                case 0
                    rem_indices = [Eindex('',1);Eindex('',1)];
                case 1
                    rem_indices = [rem_indices;Eindex('',1)];
            end
            if isempty(mute_indices);
                mute_indices = Eindex('',1);
            end
            
            reset(rem_indices);
            values = initData(T(1).data(1),rem_indices.range);
            while true
                reset(mute_indices);
                while true
                    tmp = 1;
                    for i=1:numel(T)
                        tmp = tmp*T(i).data(T(i).indices.value);
                    end
                    values(rem_indices.value) = values(rem_indices.value)+tmp;
                    if ismax(mute_indices), break; end;
                    increment(mute_indices);
                end
                if ismax(rem_indices), break; end
                increment(rem_indices);
            end
            result = Etensor(values,rem_indices);
        end
        function result = export(obj)
            result = obj.data;
        end
        function result = nabla(varargin)
            T = Etensor.empty;
            for i=1:numel(varargin)
                assert(isa(varargin{i},'Etensor'),'All inputs should be Etensors');
                T = [T;varargin{i}]; %#ok<AGROW>
            end
            [all_idx,gather_vect,scatter_vect] = unique(vertcat(T.indices),'stable');
            occurences = accumarray(scatter_vect,1,size(gather_vect));
            mute_indices = all_idx(occurences>1);
            rem_indices = all_idx(occurences==1);
            mask = ([rem_indices.range]==1);
            mute_indices = [mute_indices;rem_indices(mask)];
            rem_indices = rem_indices(~mask);
            switch(numel(rem_indices))
                case 0
                    rem_indices = [Eindex('',1);Eindex('',1)];
                case 1
                    rem_indices = [rem_indices;Eindex('',1)];
            end
            if isempty(mute_indices);
                mute_indices = Eindex('',1);
            end
            
            reset(rem_indices);
            values = initData(T(1).data(1),rem_indices.range);
            while true
                reset(mute_indices);
                while true
                    tmp = 1;
                    for i=1:(numel(T)-1)
                        tmp = tmp*T(i).data(T(i).indices.value);
                    end
                    i=i+1;
                    tmp = diff(tmp,T(i).data(T(i).indices.value));
                    values(rem_indices.value) = values(rem_indices.value)+tmp;
                    if ismax(mute_indices), break; end;
                    increment(mute_indices);
                end
                if ismax(rem_indices), break; end
                increment(rem_indices);
            end
            result = Etensor(values,rem_indices);
        end
        
        function result = plus(obj1,obj2)
            idx1 = obj1.indices.getRelevant();
            idx2 = obj2.indices.getRelevant();         
            if isempty(idx1)&&isempty(idx2)
                result = Etensor(obj1.data+obj2.data);
            else
                assert(isequal(sort(idx1),sort(idx2)),'Incompatible indices');
                result = Etensor( export(obj1) + export(obj2.permute(obj1.indices)),obj1.indices);
            end
        end
        function result = mtimes(arg1,arg2)
            result = contraction(arg1,arg2);
        end
        
        function result = voigtMatrix(obj)
            assert(isscalar(obj),'VoigtMatrix only accepts a single Etensor');
            dim = [obj.indices.range];
            assert(all(dim(1)==dim),'Voigt Notation requires that all indices have the same range');
            assert(any(dim(1)==3),'VoigtMatrix only implemented in 3D');
            
            sample0 = initData(obj.data(1),1);
            sample1 = initData(obj.data(1),1)+1;
            switch numel(obj.sz)
                case 2
                    result = repmat(sample0,6,1);
                case 4
                    result = repmat(sample0,6,6);
                otherwise
                    error('voigtMatrix implemented for tensors of order 2 or 4');
            end
            occurence = zeros(size(result));
            
            reset(obj.indices);
            while true
                [pos,factor] = obj.indices.voigtIndex(sample1);
                occurence(pos(1),pos(2)) = occurence(pos(1),pos(2))+1;
                if occurence(pos(1),pos(2))>1
                    assert(isequal(result(pos(1),pos(2)),factor*obj.data(obj.indices.value)),'Not enough symmetry to produce a Voigt notation');
                else
                    result(pos(1),pos(2)) = factor*obj.data(obj.indices.value);
                end
                if ismax(obj.indices), break; end
                increment(obj.indices);
            end
        end
        
    end
    
    methods (Static=true)
        function result = VoigtSource(data,varargin)
            sample0 = initData(data(1),1);
            sample1 = initData(data(1),1)+1;
            idx = [varargin{:}];
            assert(numel(idx)==numel(unique(idx)),'duplicate indices');
            assert(all([idx.range]==3),'voigtConstruct only implemented in 3D');
            switch numel(idx)
                case 2
                    assert(all(size(data)==[6,1]),'inconsistent data size');
                    result = Etensor(repmat(sample0,[3 3]),varargin{:});
                case 4
                    assert(all(size(data)==[6,6]),'inconsistent data size');
                    result = Etensor(repmat(sample0,[3 3 3 3]),varargin{:});
                otherwise
                    error('voigtConstruct needs a 6 by 1 matrix and 2 indices or a 6 by 6 and 4 indices');
            end
            
            reset(result.indices);
            while true
                [pos,factor] = result.indices.voigtIndex(sample1);
                result.data(result.indices.value) =  data(pos(1),pos(2))*factor^(-1);
                if ismax(result.indices), break; end
                increment(result.indices);
            end
        end
    end
end
function result = initData(sample,varargin)
if isa(sample,'SYM_BASE')
    result = SYM_ZERO(varargin{:});
    return;
end
result = zeros(varargin{:},'like',sample);
end