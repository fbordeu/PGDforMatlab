function  s = sepsqrt( arg1,Niter,res_reduc)
%SEPCOS sqrt of a separated field
%Niter = number of Newton iterations
%if negative quantities are present,... use at your own risks
%
%  s = sepsqrt( FF, [Niter=5],[res_reduc = 1e-4] )
%
%Author: Adrien Leygue.

if nargin<2
    Niter = 10;
end
if nargin<3
    res_reduc = 1e-4;
end

validateattributes(Niter,{'numeric'},{'scalar','integer','positive'},2);
validateattributes(res_reduc,{'numeric'},{'scalar','<',1,'positive'},3);

[Ndims,Ndofs,Nterms] = validateFF(arg1);

s = arg1;
r = cell(Ndims,1);
%construction du résidu = s*s-arg1
for d = 1:Ndims
    r{d} = zeros(Ndofs(d),Nterms^2+Nterms);
    for t = 1:Nterms
        r{d}(:,(1:Nterms)+(t-1)*Nterms) = bsxfun(@times,s{d},s{d}(:,t));
    end
    if d==1
        r{d}(:,Nterms^2 + (1:Nterms)) = -arg1{d};
    else
        r{d}(:,Nterms^2 + (1:Nterms)) = arg1{d};
    end
end
r = recompact(r);
err0 = sepnorm(r);

for iter = 1:Niter
    ds = sepfrac(r,s);  
    ds = sepprod(-0.5,ds);
    s = sepsum(s,ds);
    r = cell(Ndims,1);
    NtermsS = size(s{1},2);
    for d = 1:Ndims
        r{d} = zeros(Ndofs(d),NtermsS^2+Nterms);
        for t = 1:NtermsS
            r{d}(:,(1:NtermsS)+(t-1)*NtermsS) = bsxfun(@times,s{d},s{d}(:,t));
        end
        if d==1
            r{d}(:,NtermsS^2 + (1:Nterms)) = -arg1{d};
        else
            r{d}(:,NtermsS^2 + (1:Nterms)) = arg1{d};
        end
    end
    r = recompact(r);
    err = sepnorm(r);
    if (err/err0) <= res_reduc
        break;
    end
end
if (err/err0) > res_reduc
    warning('sepsqrt: maximum number of iteration reached without reaching convergence');
end

end

