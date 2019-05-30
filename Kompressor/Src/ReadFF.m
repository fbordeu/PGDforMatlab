function FF= ReadFF(filename,varargin)
%
% Function to read a semarated field from a binary file
% 
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
%
opt.single = false;

for k = 1:2:length(varargin)
    if strcmpi(varargin{k}, 'single')
        opt.single = varargin{k+1};
        continue
    end
end

FID = fopen(filename,'r');

ndims =  fread(FID,1,'int32=>int');
disp(['Number of dimension : ' num2str(ndims)])

FF = cell(ndims,1);

for i= 1:ndims
    s = fread(FID,2,'int32=>int')';
    disp(s)
    FF{i} = zeros(s);
    if opt.single
      FF{i} = reshape( fread(FID,numel(FF{i}),'float32=>double') , size(FF{i}));
    else
      FF{i} = reshape( fread(FID,numel(FF{i}),'float64=>double') , size(FF{i}));
    end
end

fclose(FID);
end