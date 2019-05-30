%Example file of the separated operations
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
%
clear all

op = 6
% 1   1+ff  
% 2   ff+ff
% 3   ff*ff
% 4   polyval(p,ff)
% 5   sin(ff)
% 6   cos(ff)
% 7   exp(ff)
% 8   max(ff)
% 9   ff^(alpha)

data = (membrane(1,50)*2*pi+2.1)/3;
[U,S,V] = svd(data,'econ');
Un = bsxfun(@times, U', diag(S))';
FF{1} = Un;
FF{2} = V;
FF = FF(:);
plot =1;

switch (op)
    case{1}
        exsol = 1+ (V*Un')';
        sepsol = sepsum(1,FF);
    case{2}
        exsol = (V*Un')' + (V*Un')';
        sepsol = sepsum(FF,FF);
    case{3}
        exsol = (V*Un')' .* (V*Un')';
        sepsol = sepprod(FF,FF);
    case{4}
        exsol =    polyval([ -0.6 -0.5 1 0],(V*Un')');
        sepsol = seppolyval([ -0.6 -0.5 1 0],FF);
    case{5}
        exsol = sin(   (V*Un')');
        sepsol = sepsin(FF);
    case{6}
        exsol = cos(   (V*Un')');
        sepsol = sepcos(FF);
    case{7}
        exsol = exp((V*Un')');
        sepsol = sepexp(FF);
    case{8}
        exsol = max(max((V*Un')));
        sepsol = sepmax(FF);
        plot = 0;
        disp(['Exact max ' num2str(exsol) ])
        disp(['Sep max ' num2str(sepsol) ])   
    case {9}
        alpha = 0.5;
        exsol = ((V*Un')').^alpha;
        tic
        sepsol = sepfracpow(FF,alpha,'around',1.3);
        toc
    otherwise
end

if plot
    figure(1)
    subplot(1,3,1)
    surf(exsol);
    title ('Exact')
    shading flat;
    
    subplot(1,3,3)
    SS = (sepsol{2}*sepsol{1}')';
    surf(SS);
    title (['Sep  with ' num2str(size(sepsol{1},2)) ' terms'])
    shading flat;
    
    subplot(1,3,2)
    surf(exsol - SS);
    title ('Error')
    shading flat;
else
    figure(1)
    clf
    surf(data);
    shading flat;
end