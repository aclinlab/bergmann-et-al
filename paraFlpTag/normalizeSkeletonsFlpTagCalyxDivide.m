function skeletons = normalizeSkeletonsFlpTagCalyxDivide(filenames, spacing)

% Let's start with filenames as an cell array of strings specifying the .mat
% files

% calling function must supply the value of 'spacing' which is the node
% spacing on the 'standardized' skeleton

% The purpose of this function is normalize the lengths of the vertical
% lobe, peduncle, and horizontal lobes to single "standard" lengths
% each file will be a .mat file containing multiple activityMap objects -
% they would all have the same skeleton but would be different movies (e.g.
% odor, odor+ATP, ATP alone, and the pipette in different locations)

% Let's try returning the array of skeleton objects ('skeletons')

numFiles = length(filenames);

% length will store the length of each branch in each file
lengths = zeros(numFiles,4); % 1 = c, 2 = p, 3 = h, 4 = v
% the length
peduncleLength = zeros(numFiles);
% junctionIndices will store the index of the node on each skeleton that is
% the junction
junctionIndices = zeros(numFiles,1);

for i=1:numFiles
    i
    % loop through all files
    %  read in .mat file
    load(filenames{i}); % assumes the .mat file contains a variable called s which has the skeleton
    skeletons(i) = s;
    skeletons(i) = skeletons(i).divideSkeleton; % use the normal vector to divide the calyx from the peduncle
    
    %   find the node that has >2 connections - that's the junction
    nodeLinks = {skeletons(i).nodes.links};
    numLinksPerNode = zeros(length(nodeLinks),1);
    for j=1:length(nodeLinks)
        numLinksPerNode(j) = length(nodeLinks{j});
    end
    
    [maxLinks,junctionIndices(i)] = max(numLinksPerNode);
    if maxLinks~=3
        error('the node with the most links does not have 3 links!');
    end
    % to measure its peduncle, vertical lobe and
    %  horizontal lobe
    %   just measure the distance along the skeleton between each end point and
    %   the junction
    % Whereas in Amin et al. 2020, the skeleton is always drawn from vertical lobe tip to
    % calyx, then the horizontal lobe is added from junction to horizontal
    % lobe tip last [so the endNodes were in the order: vertical, calyx,
    % horizontal]
    % In contrast, here we have drawn the skeletons first calyx to
    % horizontal lobe, then the vertical lobe is added from the junction to
    % the vertical tip.
    % Thus endNodes are in the order: calyx, horizontal, vertical
    endNodes = skeletons(i).findEndNodes();

    calyx = endNodes(1);
    horEnd = endNodes(2);
    vertEnd = endNodes(3);
    % Get the distances from the calyx end to all nodes in the skeleton
    [distancesFromCalyx,~] = skeletons(i).getDistances(calyx);

    % Get the distances from the junction node to all nodes in the skeleton
    [distancesFromJunction,~] = skeletons(i).getDistances(junctionIndices(i));
    
    if length(endNodes)~=3
        error('number of endNodes ~= 3!');
    end

    calyxLength = distancesFromCalyx(skeletons(i).dividingNode); % the dividingNode-th node is the division between calyx and peduncle
    peduncleLength = distancesFromJunction(skeletons(i).dividingNode);
    horLength = distancesFromJunction(horEnd);
    vertLength = distancesFromJunction(vertEnd);
    lengths(i,:) = [calyxLength, peduncleLength, horLength, vertLength];
    
end
meanLengths = mean(lengths,1);
disp('mean lengths =');
meanLengths

% Divide every length by the mean length for that lobe
stretchFactors = lengths./repmat(meanLengths, numFiles, 1);

% Loop over all files again
% For each set of 9 skeletons (each file), modify the skeleton
for i=1:numFiles
    % label the branches (divide the skeleton into segments)
    % here it should be calyx = 1, peduncle = 2, horizontal = 3, vertical = 4
    % Note this is different from Amin et al 2020 or from the other
    % normalize function which doesn't divide calyx from peduncle
    skeletons(i) = skeletons(i).labelBranchesSepCalyx();
    skeletons(i).drawDividedSkeleton; %%%%%%%%%%%%%%%%%%%COMMENT OUT - FOR DEBUGGING ONLY
    
    % Re-draw the  spaced nodes, emanating from the junction
    % Set the starting point to be the junction
    skeletons(i).userStartingPoint = skeletons(i).nodes(junctionIndices(i)).realCoords;
    % skeletons(i).userStartingPoint = skeletons(i).nodes(1).realCoords; % for debugging
    skeletons(i) = skeletons(i).createSpacedNodes(spacing*stretchFactors(i,:));
    skeletons(i) = skeletons(i).createVoronoiMask();
end

end