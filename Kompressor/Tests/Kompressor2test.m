function res = Kompressor2test( input_args )
%test kompressor with the new XdmfGrid
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
%
res = true;

%data1 = readpxdmf('data/ParaMat3D.pxdmf');
data2 = readpxdmf2([ 'data' filesep() 'ParaMat3D.pxdmf'] );
data2R = recompact(data2,'verbose',true);
end

