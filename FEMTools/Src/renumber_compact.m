function [result, N] = renumber_compact(arg)
    %replaces the elements of arg by their numbering should they be sorted
    %duplicate elements receive the same number
    %N is the maximum of result.
    %
    % This file is subject to the terms and conditions defined in
    % file 'LICENSE.txt', which is part of this source code package.
    %
    % Principal developer : Adrien Leygue (Adrien.Leygue@ec-nantes.fr)
    %
    [~ ,gather_vect,scatter_vect] = unique(arg,'sorted');
    result = reshape(scatter_vect,size(arg));
    N = numel(gather_vect);
end

