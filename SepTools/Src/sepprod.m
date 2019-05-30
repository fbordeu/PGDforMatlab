function  s = sepprod( varargin )
%SEPPROD product of 2 or more separated field
%   
% s = sepprod( FF1, FF2, ..., FFn, options, ... )
%
% FFi can be a scalar or a separated field
%
% See also PushRecompactOptions, PopRecompactOptions.
% 
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
%

% if allNterms is bigger than dNterms a recursion is used to keep memory
% foot print low.

  [args1, args2] = SplitInputs(varargin);
  
  s = args1{1};
  for i = 2:numel(args1)
    s = sepprod_internal(s,args1{i},args2);
  end
  
end


function  s = sepprod_internal( arg1,arg2,recompactoptions )

% scalar multiplication
if ~iscell(arg1)
    s = arg2;
    s{1} = arg2{1}*arg1;
    return 
end

if ~iscell(arg2)
    s = arg1;
    s{1} = arg1{1}*arg2;
    return 
end

[Ndims,Ndofs,allNterms] = multiValidateFF(arg1,arg2);   

% field multiplication
s = cell(Ndims,1);

nmodes1 = allNterms(1);
nmodes2 = allNterms(2);

% 2 G divided by the size of one terms
dNterms = round(sqrt(20.e9/(8*Ndims*sum(Ndofs))));

if all(allNterms<=dNterms)

elseif allNterms(1)>dNterms
    limit = floor(allNterms(1)/2);
    tmp = cellfun(@(x) x(:,1:limit),arg1,'uniformoutput',false);
    s = sepprod_internal(tmp,arg2,recompactoptions);
    tmp = cellfun(@(x) x(:,(limit+1):allNterms(1)),arg1,'uniformoutput',false);
    s = sepsum(s,sepprod_internal(tmp,arg2,recompactoptions),recompactoptions{:});
    return
elseif allNterms(2)>dNterms
    limit = floor(allNterms(2)/2);
    tmp = cellfun(@(x) x(:,1:limit),arg2,'uniformoutput',false);
    s = sepprod_internal(arg1,tmp,recompactoptions);
    tmp = cellfun(@(x) x(:,(limit+1):allNterms(2)),arg2,'uniformoutput',false);
    s = sepsum(s, sepprod_internal( arg1,tmp,recompactoptions),recompactoptions{:});
    return
end
     
for d = 1:Ndims
  size1 = size(arg1{d},1);
  mat = zeros(size1,nmodes1*nmodes2);
    cpt =1;
    for i = 1: nmodes1
      for j = 1:nmodes2
        mat(:,cpt) = arg1{d}(:,i).*arg2{d}(:,j);
        cpt = cpt+1;
      end
    end
    s{d} = mat;
end

s = recompact(s,recompactoptions{:});



end

