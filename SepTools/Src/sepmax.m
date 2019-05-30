function [Y,I] = sepmax( arg1, varargin )
%SEPMAX max of a separated field
%
% [Y,I] = sepmax( FF, options )
%
%  NOTE: The field MUST be positive 
%   FF always > 0  
%
%  options are options for the internal recompact.
%  and/or 
%  Number of max searched in the field ('Nmax',1)
%  Maximal number of iteration for the power algo ('Mmax_iters',100)
%
% Y the maximus
% I the index
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
%
[Ndims,Ndofs,Nterms] = validateFF(arg1);

Nmax =1;
Mmax_iters = 100;

for k = 1:length(varargin)
    if strcmpi(varargin{k}, 'Nmax')
        Nmax = varargin{k+1};
        varargin = varargin([1:(k-1) (k+2):end]);
        break
    end
end

for k = 1:length(varargin)
    if strcmpi(varargin{k}, 'Mmax_iters')
        Mmax_iters = varargin{k+1};
        varargin = varargin([1:(k-1) (k+2):end]);
        break
    end
end

Y = zeros(Nmax,1);
I = zeros(Nmax,Ndims);

% for every max point
cpt =1;
while (cpt < Nmax+1 ) 
    
    prod = internal_sepmax( arg1, Nmax, Mmax_iters, varargin{:});
    
    Is = zeros(1,Ndims);
        
    for i = 1:Ndims
        [~,Is(i)] = max(abs(prod{i}(:,1)));
    end

    Y(cpt) = sepeval(arg1,Is);
    I(cpt,:)  =  Is;
    
    fxn = zeros(Ndofs(1),1);
    %disp(['V : ' num2str(Y(cpt)) ' P : '  num2str(Is)])
    fxn(Is(1),1) = -1*Y(cpt);
    %keyboard
    arg1{1} = [arg1{1} fxn ];
    for i = 2:Ndims
    	fxn = zeros(Ndofs(i),1);
        fxn(Is(i),1) = 1;
        arg1{i} = [arg1{i} fxn ];
    end
    
    %figure(2)
    %clf
    %SS = (arg1{1}*arg1{2}')';
    %surf(SS); 
    %shading flat
    %title(num2str(cpt))
    %hold on
    %plot3(I(1:(cpt),1),I(1:(cpt),2),Y(1:(cpt)),'x');
    %pause(0.01)
    
    cpt = cpt +1;
    
end

end



function prod = internal_sepmax( arg1,Nmax,Mmax_iters,varargin )
dims1 =   numel(arg1);
prod = arg1;
cpt = 0;
while(true)
    
    % product
    %  prod = prod*arg1 <- this is robust but slow
    %  prod = prod*prod <-  this is really fast but not robust
    %  prod = prod*(prod+arg1)  < - this is a mix
    prod = sepprod(prod,sepsum(prod,arg1,'recompact',false), 'fp_tol',1e-5,'res_reduc',1e-4,'verbose',0,varargin{:});

    
    %figure(3)
    %normalization
    for i = 1:dims1
      alphas = sqrt(sum(prod{i}.*prod{i}));
      if(i == 1); alphas1 = alphas ; end;
      %if (i==1) semilogy(alphas);end
      %hold on
      nor = 1/norm(alphas);  
      prod{i} = prod{i}*nor; 
    end
    
    
    %if mod(cpt,1) == 0
    %    figure(1)
    %    SS = (prod{2}*prod{1}')';
    %    surf(SS);
    %    shading flat
    %    pause(0.01)
    %end
    
    %disp(size(prod{1}))
    cpt = cpt +1;
 
    %figure(3)
    %plot(prod{1})
    %disp(alphas1)
    %cpt
    %semilogy(alphas1)
    %if( alphas1(end-1) > 1e5*alphas1(end) ); break; end;
    if (size(prod{1},2) <= 3 || cpt > Mmax_iters) ; break;end
    if (size(prod{1},2) == 1 )
       a = 1;
       for i = 1:dims1   
          [~,index] = max(abs(prod{index}));
          a = a*max(prod{i}(i)); 
       end
       if a > 0.5; break; end;
    end
end

end


