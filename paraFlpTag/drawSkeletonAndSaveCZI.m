nodeSpacing = 10; % microns

% load mask data
% [maskFileName, maskPathName] = uigetfile('*.mat', 'Pick the .mat file with your skeleton object');
% load(strcat(maskPathName, filesep, maskFileName),'mask','bkgndmask','fileName','pathName');
% 
% %% read file
% 
% % Use this line if your 2 channels are in one file
[stack, rowCal, columnCal, sliceCal, fileName, pathName] = readCZI2channels();
% % Use this line if your 2 channels are split across 2 files
% [stack, rowCal, columnCal, sliceCal, fileName, pathName] = readCZI2channels2files();
% 

%% user draws the skeleton endpoints
disp('Click on the window to define the skeleton nodes');
disp('Rotate by pressing j, k, i or m');
disp('Click twice: first to define a line in 3d space, then rotate the picture,');
disp(' then click where along the line defined by the first click you want the point to be');
disp('After you defined one node, press n for next to define the next node');
disp('Or press t for terminate to say you''re finished');
s = skeleton(mask,[rowCal columnCal sliceCal]);

%% user defines starting point
disp('Now click on the window to define the starting point of the skeleton');
s = s.defineStartingPoint(mask);

timeStamp = datestr(now,'yymmddTHHMMSS');
% save in case of crash
save(strcat(fileName,'-skeleton',timeStamp,'.mat'),'s','bkgndmask');

% %% this requires no user input
s = s.createSpacedNodes(nodeSpacing);
s = s.createVoronoiMask;
s = s.trimVoronoiMaskEndpoints; % trim the calyx off of the voronoi masks

channel1 = stack(:,:,:,1);
channel2 = stack(:,:,:,2);

% find if there is any background mask in each z-slice
isBkgndMaskInThisSlice = max(bkgndmask,[],[1 2]);

for i=1:size(channel1,3)
    channel1ThisFrame = channel1(:,:,i);
    channel2ThisFrame = channel2(:,:,i);
    if isBkgndMaskInThisSlice(i)
        channel1(:,:,i) = channel1ThisFrame - mean(channel1ThisFrame(bkgndmask(:,:,i)),'all');
        channel2(:,:,i) = channel2ThisFrame - mean(channel2ThisFrame(bkgndmask(:,:,i)),'all');
    else
        % find the nearest slice where there is bkgnd
        existingBkgndSlices = find(isBkgndMaskInThisSlice);
        [~,distToNearestbkgnd_idx] = min(abs(i-existingBkgndSlices));
        nearestBkgndSlice = existingBkgndSlices(distToNearestbkgnd_idx);
        channel1(:,:,i) = channel1ThisFrame - mean(channel1ThisFrame(bkgndmask(:,:,nearestBkgndSlice)),'all');
        channel2(:,:,i) = channel2ThisFrame - mean(channel2ThisFrame(bkgndmask(:,:,nearestBkgndSlice)),'all');
    end
end

%% this would be a bad way to measure background if there is other bright stuff in the background
%bkgnd1 = mean(channel1(~mask),'all');
%bkgnd2 = mean(channel2(~mask),'all');

%% calculate intensity at each voronoi volume
% bkgnd1 = mean(channel1(bkgndmask),'all');
% bkgnd2 = mean(channel2(bkgndmask),'all');
numNodes = length(s.spacedNodes);
[c1, c2] = deal(zeros(numNodes,1));
for i=1:numNodes
    c1(i) = mean(channel1(s.voronoiMask(:,:,:,i)),'all'); % - bkgnd1;
    c2(i) = mean(channel2(s.voronoiMask(:,:,:,i)),'all'); % - bkgnd2;
end

%% Check if the skeleton is OK
s.drawSkeleton(c1./c2)

%save in case of crash
save(strcat(fileName,'-skeleton',timeStamp,'.mat'),'s','c1','c2','bkgndmask');

%% calculate intensity at each half
s = s.divideMask;
[c1split, c2split] = deal(zeros(2,1));
for i=1:2
    c1split(i) = mean(channel1(s.splitMask(:,:,:,i)),'all'); % - bkgnd1;
    c2split(i) = mean(channel2(s.splitMask(:,:,:,i)),'all'); % - bkgnd2;
end

%% check if split in half is OK
s.drawSplitSkeleton(c1split./c2split);

%% save data
save(strcat(fileName,'-skeleton',timeStamp,'.mat'),'s','c1','c2','c1split','c2split','bkgndmask')
