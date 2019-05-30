function  s = sepsum( varargin )
%SEPSUM Sum of two or more separated fields with internal recompact
%
% function  s = sepsum( FF1, FF2,..., FFn, options, ... )
%
% FFi can be a scalar or a field 
%
%  options are options for the internal recompact.
%
%  Note if a scalar is used, no recopaction is applied
%
% See also PushRecompactOptions, PopRecompactOptions, recompact.
%
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
%
  [args1, args2] = SplitInputs(varargin);
  
  s = args1{1};
  for i = 2:numel(args1)
    s = sepsum_internal(s,args1{i},args2);
  end
  
end

function  s = sepsum_internal( arg1,arg2,recompactoptions)

if iscell(arg1)
  if  iscell(arg2)
    % Validation of the inputs
    multiValidateFF(arg1,arg2);   
    
    % addition
    s = cellfun(@(a1,a2)([a1 a2]),arg1,arg2,'UniformOutput',false);
    
    % recompaction
    s = recompact(s,recompactoptions{:});
    
  else
    % we treat arg2 as a scalar
    s = sepsum_internal(arg2,arg1,recompactoptions );
  end
else
  if  iscell(arg2)  
    % we treat arg1 as a scalar
    
    if arg1 == 0 
        s = arg2;
        return
    end
    
    dims2 =   numel(arg2);
    s = cell(dims2,1);
  
    for d = 1:dims2
        s{d} = ones(size(arg2{d},1),1 );
    end
  
    s{1} = s{1}*arg1;
    s = sepsum_internal( s,arg2,{recompactoptions{:} 'recompact' false}  );
  else
      
    % in the case arg1 and arg2 are scalars
    s = arg1+arg2;
  end
end

end


    


