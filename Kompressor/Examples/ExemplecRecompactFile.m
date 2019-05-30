%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
%

close all
clear all

% Generator of the data
load('durer','X','map');
[U,S,V] = svd(double(X),'econ');
Un = bsxfun(@times,U,diag(S)');


% put the date in the structure
FF{1} = Un;
FF{2} = V;

% Write data to the disk
DumpFF('Data.bin',FF);
DumpFF('Data_single.bin',FF,'single', true);


% recompact the data

%system('gdb -ex r -args ../Src/cRecompactFile -MAM 20 Data.bin')
system('../Src/cRecompactFile_Linux -MAM 20 Data.bin');
%system('gdb -ex r -args ../Src/cRecompactFile_Darwin -MAM 20 +float Data_single.bin')
system('../Src/cRecompactFile_Linux -MAM 20 -lastImproveModesLoop 10 +float Data_single.bin');

% Verification of the data from the disk 
disp('Diff between  write and read')
FF_save = ReadFF('Data.bin');
disp([ norm(FF_save{1} - FF{1}) norm(FF_save{2} - FF{2}) ] )


FF_new = ReadFF('out_Data.bin');
FF_new_single = ReadFF('out_Data_single.bin','single',true);


figure (1);
subplot(2,3,2)
imagesc((V*Un')')
colormap(map)
axis equal
title 'original'

subplot(2,3,1)
n =(FF_new{2}*FF_new{1}')';
imagesc(n)
colormap(map)
axis equal
title 'double 20 terms'


subplot(2,3,3)
imagesc((FF_new_single{2}*FF_new_single{1}')')
colormap(map)
axis equal
title 'single 20 terms'

subplot(2,3,5)
dif = (FF_new{2}*FF_new{1}')' - (FF_new_single{2}*FF_new_single{1}')';
imagesc(dif)
colormap(map)
axis equal
title(['dif ' num2str(max(dif(:))/max(n(:)) )])


