function DumpFF(filename,FF,varargin)
%
% Function to write a serparated field from a binary file
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

FID = fopen(filename,'w');

ndims = numel(FF);
disp(['Number of dimension : ' num2str(ndims)])

cout = 0;
cout = cout + fwrite(FID,int32(ndims),'int32');
for i= 1:ndims
    cout = cout + fwrite(FID,int32(size(FF{i})),'int32');
    if opt.single 
        cout = cout + fwrite(FID,double(FF{i}),'float');
    else
        cout = cout + fwrite(FID,double(FF{i}),'double');
    end
end

%for i= 1:ndims
%end

fclose(FID);