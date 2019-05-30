function result =conv2SYM(arg)
    %Converts numeric argument to their symbolic counterpart.
    %
    % This file is subject to the terms and conditions defined in
    % file 'LICENSE.txt', which is part of this source code package.
    %
    % Principal developer : Adrien Leygue (Adrien.Leygue@ec-nantes.fr)
    %
    if isa(arg,'SYM_BASE')
        result = arg;
    else
        result = SYM_REAL(arg);
        for i=1:numel(result)
            switch result(i).value
                case 0
                    result(i) = SYM_ZERO;
                case 1
                    result(i) = SYM_ONE;
            end
        end
    end
end