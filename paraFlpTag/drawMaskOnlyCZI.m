%% settings
sliceInterval = 3;

scaleIntensity = 1;

%% read file
[stack, rowCal, columnCal, sliceCal, fileName, pathName] = readCZI2channels();

%% user draws the mask
disp(fileName);
disp('Outline the mushroom body in each slice');
disp(strcat('It is only displaying one slice out of every ', num2str(sliceInterval),' slices'));
disp('The software will interpolate the mask between the intervening slices');
mask = skeleton.outlineObject(stack,2,sliceInterval,scaleIntensity);

%% measure background using hand drawn mask
disp('Draw an outline to define the ROI of the background fluorescence');
disp('use somewhere that is black for both the green and red channel');
disp('ideally somewhere deep in the stack');
disp('You can draw ROIs in more than one slice');
disp('Because the software will interpolate between slices, dont'' draw a single');
disp('ROI to go across multiple slices');
bkgndmask = skeleton.outlineObject(stack,1,sliceInterval,scaleIntensity);
timeStamp = datestr(now,'yymmddTHHMMSS');
save(strcat(fileName,'-maskData',timeStamp,'.mat'),'mask','bkgndmask','fileName','pathName');

