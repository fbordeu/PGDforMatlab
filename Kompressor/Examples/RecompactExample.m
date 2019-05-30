%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
%

close all
%clear all

% Read data 
street1=imread('street1.jpg');
street1=double(street1);

% SVD 
[U,S,V] = svd((street1(:,:,1)+street1(:,:,2)+street1(:,:,3))/10,'econ');

% alphas
Un = bsxfun(@times, U', diag(S))';

%figure
figure(1)
subplot(1,2,1)
image((V*Un')');
colormap(gray)


% we put the data in the structure
FF{1} = Un;
FF{2} = V;

% Recompact
data = recompact(FF,'max_added_modes',10,'verbose',1,'lastImproveModesLoop',20);

% for this you must launch the cRecompactServer fisrt
data1 = recompact(FF,'UseTCPIP',true,'host','localhost','port','5556','max_added_modes',20,'verbose',1,'lastImproveModesLoop',20);
data20 = recompact(FF,'UseTCPIP',true,'host','localhost','port','5556','max_added_modes',20,'verbose',1,'lastImproveModesLoop',1);

%figure
figure(1)
subplot(1,2,2)
image((data{2}*data{1}')');
colormap(gray)
pause(0.1)



