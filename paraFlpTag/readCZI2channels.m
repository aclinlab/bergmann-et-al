function [stack, rowCal, columnCal, sliceCal, fileName, pathName] = readCZI2channels(varargin)

% Assumes that there are 2 channels

if nargin==2
    fileName = varargin{1};
    pathName = varargin{2};
else
    [fileName, pathName] = uigetfile('*', 'Pick a file');
end

% Get the extension
[~,~,extension] = fileparts(fileName);

disp('Now reading in the file');

data = bfopen(strcat(pathName,fileName));
numImages = size(data{1},1);
numSlices = numImages/2;
XYsize = size(data{1}{1,1});
stack = zeros(XYsize(1), XYsize(2), numSlices, 2);
for i=1:numSlices
    stack(:,:,i,1) = data{1}{(i-1)*2+1,1};
    stack(:,:,i,2) = data{1}{i*2,1};
end

% swap the order of the channels for .czi files
% this is a hack for Katie's data where the channels are backwards on .czi
% files
if strcmp(extension,'.czi')
    stack = stack(:,:,:,[2 1]);
end

omeMeta = data{1,4};
rowCal = double(omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER));
columnCal = double(omeMeta.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROMETER));
sliceCal = double(omeMeta.getPixelsPhysicalSizeZ(0).value(ome.units.UNITS.MICROMETER));

end
