[textFileName, textFilePathName] = uigetfile('*.txt', 'Pick txt file specifying input');
fid = fopen(strcat(textFilePathName,textFileName));
% filename pathname conditions 
textFileInput = textscan(fid, '%s %s %s %s %s', 'Delimiter', '\t', 'HeaderLines', 1);
fclose(fid);

numFiles = length(textFileInput{1});
fullPathFileNames = cell(numFiles,1);
for i=1:numFiles
    thisFileName = textFileInput{1}{i};
    % remove double quotes that Matlab inserts in file names with special
    % characters, like commas
    if thisFileName(1)=='"'
        thisFileName = thisFileName(2:(end-1));
    end
    % concatenate pathname with filename
    fullPathFileNames{i} = strcat(textFileInput{2}{i}, thisFileName);
end
fileNames = textFileInput{1};
pathNames = textFileInput{2};
conditions = textFileInput{3};
movieFileNames = textFileInput{4};
moviePathNames = textFileInput{5};

%filenames is now a cell array with the full path of each .mat file

spacing = 10;
% normalizeSkeletons will read in each .mat file to get its skeleton, then
% stretch or squish the spaced nodes so they will all line up together
normalizedSkeletons = normalizeSkeletonsFlpTagCalyxDivide(fullPathFileNames, spacing);
disp('trimming endpoints')
for p=1:length(normalizedSkeletons)
    % Cludge - we want the code to assign the junction to be on the same branch
    % for all mat files, but for some mat files it assigns the junction to
    % different branches. This sets the first node to be on branch 2, i.e. peduncle.
    normalizedSkeletons(p).spacedNodes(1).branchNumber=2;
    % trim off the ends of the mask that fall outside the manually drawn
    % skeleton
    normalizedSkeletons(p) = normalizedSkeletons(p).trimVoronoiMaskEndpoints;
end

[branches, branches_c1, branches_c2] = deal(cell(numFiles,4)); %4 branches per skeleton
[mergeBranches, mergeBranches_c1, mergeBranches_c2] = deal(cell(2, numFiles));
results = zeros(numFiles,length(normalizedSkeletons(1).spacedNodes));

for i=1:numFiles
%     % remove the tag '-skeleton.mat' to recover the original name of the
%     % .tiff file
%     % the .tiff file must be in the same director as the '-skeleton
%     index = strfind(fileNames{i},'-skeleton.mat');
%     movieFileName = fileNames{i}(1:(index-1));

    disp('Now reading in the  file');

    %% read file
    [stack, rowCal, columnCal, sliceCal, fileName, pathName] = readCZI2channels(movieFileNames{i},moviePathNames{i});

    channel1 = stack(:,:,:,1);
    channel2 = stack(:,:,:,2);

    % calculate normalised GFP/Cherry at each voronoi volume
    load(fullPathFileNames{i},'bkgndmask'); % load the bkgndmask from the -skeleton.mat file

    % find if there is any background mask in each z-slice
    isBkgndMaskInThisSlice = max(bkgndmask,[],[1 2]);

    for j=1:size(channel1,3)
        channel1ThisFrame = channel1(:,:,j);
        channel2ThisFrame = channel2(:,:,j);
        if isBkgndMaskInThisSlice(j)
            channel1(:,:,j) = channel1ThisFrame - mean(channel1ThisFrame(bkgndmask(:,:,j)),'all');
            channel2(:,:,j) = channel2ThisFrame - mean(channel2ThisFrame(bkgndmask(:,:,j)),'all');
        else
            % find the nearest slice where there is bkgnd
            existingBkgndSlices = find(isBkgndMaskInThisSlice);
            [~,distToNearestbkgnd_idx] = min(abs(j-existingBkgndSlices));
            nearestBkgndSlice = existingBkgndSlices(distToNearestbkgnd_idx);
            channel1(:,:,j) = channel1ThisFrame - mean(channel1ThisFrame(bkgndmask(:,:,nearestBkgndSlice)),'all');
            channel2(:,:,j) = channel2ThisFrame - mean(channel2ThisFrame(bkgndmask(:,:,nearestBkgndSlice)),'all');
        end
    end

    %bkgnd1 = mean(channel1(bkgndmask),'all');
    %bkgnd2 = mean(channel2(bkgndmask),'all');
    numNodes = length(normalizedSkeletons(i).spacedNodes);
    [c1, c2] = deal(zeros(numNodes,1));
    for j=1:numNodes
        c1(j) = mean(channel1(normalizedSkeletons(i).voronoiMask(:,:,:,j)),'all');% - bkgnd1;
        c2(j) = mean(channel2(normalizedSkeletons(i).voronoiMask(:,:,:,j)),'all');% - bkgnd2;
    end
    result = c1./c2;
    results(i,:) = result;

    % calyx = 1, peduncle = 2, horizontal = 3, vertical = 4
    branchNumbers = [normalizedSkeletons(i).spacedNodes.branchNumber];
    for k=1:size(branches,2)
        branches{i,k} = squeeze(result(branchNumbers==k));
        branches_c1{i,k} = squeeze(c1(branchNumbers==k));
        branches_c2{i,k} = squeeze(c2(branchNumbers==k));
    end

    % 1: calyx-->vertical (branch 1 inverted, then branch 2 inverted, then 4)
    mergeBranches{1,i} = [flipud(branches{i,1}); flipud(branches{i,2}); branches{i,4}];
    mergeBranches_c1{1,i} = [flipud(branches_c1{i,1}); flipud(branches_c1{i,2}); branches_c1{i,4}];
    mergeBranches_c2{1,i} = [flipud(branches_c2{i,1}); flipud(branches_c2{i,2}); branches_c2{i,4}];
    % 2: calyx-->horizontal (branch 1 inverted, then branch 2 inverted, then 3)
    mergeBranches{2,i} = [flipud(branches{i,1}); flipud(branches{i,2}); branches{i,3}];
    mergeBranches_c1{2,i} = [flipud(branches_c1{i,1}); flipud(branches_c1{i,2}); branches_c1{i,3}];
    mergeBranches_c2{2,i} = [flipud(branches_c2{i,1}); flipud(branches_c2{i,2}); branches_c2{i,3}];



end

branchKey = {'CtoV','CtoH'};
for j=1:length(branchKey)
    csvwrite([branchKey{j}, '.csv'], cell2mat(mergeBranches(j,:)));
    csvwrite([branchKey{j}, '_c1.csv'], cell2mat(mergeBranches_c1(j,:)));
    csvwrite([branchKey{j}, '_c2.csv'], cell2mat(mergeBranches_c2(j,:)));
end

timeStamp = datestr(now,'yymmddTHHMMSS');
save(strcat('normalizedSkeletons-',timeStamp,'.mat'), 'normalizedSkeletons', 'branches');