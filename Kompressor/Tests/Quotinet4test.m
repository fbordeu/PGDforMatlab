function res = quotient4test
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
%
clear all
n = 100;
% 
[U,S,V]  = svd((membrane(1,50)));
U = U(:,1:n);
S = S(1:n,1:n);
V = V(:,1:n);

% alphas
Un = bsxfun(@times, U', diag(S))';

FF1 =  {Un; V};

n2 = n;
[U,S,V]  = svd(membrane(1,50)+2);
U = U(:,1:n2);
S = S(1:n2,1:n2);
V = V(:,1:n2);

% alphas
Un = bsxfun(@times, U', diag(S))';

FF2 = {Un; V};

%quotient({},{},'recompile',true)

tic
solv6   = quotient(FF1,FF2,'max_added_modes',100,'useV6',true,'fp_max_iter',40,'res_reduc',1e-6,'improve_modes',true);
v6toc = toc;

tic
solEasy = quotient(FF1,FF2,'max_added_modes',100,'useEasyQuotient',true,'fp_max_iter',40,'res_reduc',1e-6);
EasyQuotienttoc = toc;

tic
solc    = quotient(FF1,FF2,'max_added_modes',100,'useV6',false,'fp_max_iter',40,'res_reduc',1e-6,'improve_modes',true);
ctoc = toc;
%sol1 =recompact (FF1);
[ v6toc EasyQuotienttoc ctoc ]

clf
figure(1);
subplot(3,3,1)
surface((FF1{1}*FF1{2}')./(FF2{1}*FF2{2}'))
title('Exact')

subplot(3,3,3)
surface(solv6{1}*solv6{2}')
title([ 'SolV6 (' num2str(size(solv6{1},2)) '): ' num2str(v6toc)  ])

subplot(3,3,2)
surface((FF1{1}*FF1{2}')./(FF2{1}*FF2{2}') - solv6{1}*solv6{2}')
res = max(max((FF1{1}*FF1{2}')./(FF2{1}*FF2{2}') - solv6{1}*solv6{2}'));

subplot(3,3,7)
surface(solEasy{1}*solEasy{2}')
title([ 'solEasy (' num2str(size(solEasy{1},2)) '): ' num2str(EasyQuotienttoc) ])

subplot(3,3,4)
surface((FF1{1}*FF1{2}')./(FF2{1}*FF2{2}') - solEasy{1}*solEasy{2}')
res = max(max(max((FF1{1}*FF1{2}')./(FF2{1}*FF2{2}') - solEasy{1}*solEasy{2}')),res);

subplot(3,3,9)
surface(solc{1}*solc{2}')
title([ 'solc (' num2str(size(solc{1},2)) '): ' num2str(ctoc) ])

subplot(3,3,5)
surface((FF1{1}*FF1{2}')./(FF2{1}*FF2{2}') - solc{1}*solc{2}')

res = max(max(max((FF1{1}*FF1{2}')./(FF2{1}*FF2{2}') - solc{1}*solc{2}')),res);

res = res< 1e-5;
