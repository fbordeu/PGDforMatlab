function res = SepToolstests()
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
%
disp('Running SepToolstests')
res = true;

data = (membrane(1,50)*2*pi+3)/3;
[U,S,V] = svd(data,'econ');
Un = bsxfun(@times, U', diag(S))';
FF= { Un V}';

F = (V*Un');

sex = 1 + F;
sepsol = sepsum(1,FF);
res = res & MaxError('sum',sepsol, sex);

sex = F + F;
sepsol = sepsum(FF,FF);
res = res & MaxError('sum',sepsol, sex);

sex = F .* F;
sepsol = sepprod(FF,FF);
res = res & MaxError('prod',sepsol, sex);

sex = polyval([ -0.6 -0.5 1 0],F);
sepsol = seppolyval([ -0.6 -0.5 1 0],FF);
res = res & MaxError('polyval',sepsol, sex);

sex = sin(F);
sepsol = sepsin(FF);
res = res & MaxError('sin',sepsol, sex);

sex = cos(F);
sepsol = sepcos(FF);
res = res & MaxError('cos',sepsol, sex);


sex = exp(F);
sepsol = sepexp(FF);
res = res & MaxError('exp',sepsol, sex);

sex = max(F(:));
sepsol = sepmax(FF);
res = res & MaxError('max',sepsol, sex);

alpha = 0.5;
sex = (F).^alpha;
sepsol = sepfracpow(FF,alpha,'around',2);
res = res & MaxError('pow',sepsol, sex);

end

function res = MaxError(name,sepsol, sex)

if isa(sepsol,'numeric') && isa(sex,'numeric')
    err = abs(sepsol-sex);
else
    err = max(max(abs(sepsol{2}*sepsol{1}' - sex)));    
end


if  err > 1e-4 
    disp([ 'Error in ' name ])
    disp([ 'Error more than ' num2str(1e-4) ])
    res = false;
else
    res = true;
end

    
end