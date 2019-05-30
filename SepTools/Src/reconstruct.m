function result = reconstruct(FF)
%
% Take a FF and reconstruct the field in full
%
%
    [Ndims,Ndofs,~] = validateFF(FF);
    
    P1 = my_khatriraoR(FF(1:floor(Ndims/2)));
    P2 = my_khatriraoR(FF((floor(Ndims/2)+1):end));
    result = P1*P2';
    result = reshape(result(:),Ndofs');
    
end

function P = my_khatriraoR(FF)
tmp = numel(FF);
    switch tmp
        case 1
            P = FF{1};
        case 2 
            P = khatrirao(FF,'r');
        otherwise
            P1 = my_khatriraoR(FF(1:floor(tmp/2)));
            P2 = my_khatriraoR(FF((floor(tmp/2)+1):end));
            P = khatrirao({P1,P2},'r');
    end
end