function ROIs = clickCircles(figHandle, axHandle, circleSize, varargin)
% clickCircles:  Interactively specify multiple circles
% varargin is if you want to give existing ROIs
% note the existing mask must be 3d (x, y, # ROIs)!
%
% Overlays an imfreehand ROI on an image.
% Gives the ability to tweak the ROI by adding and
%   subtracting regions as needed, while updating
%   the ROI boundaries as an overlay on the image.
% Returns a logical matrix of the same size as the
%   overlain image.
%
% Requires alphamask:
%   http://www.mathworks.com/matlabcentral/fileexchange/34936
%
% Usage:
%   bwMask = fhroi(figHandle, [axHandle])
%     axHandle: handle to axes on which to operate (optional)
%       bwMask: ROI mask as logical matrix
%
% Example:
%   f = figure;
%   I = rand(20) + eye(20);
%   imshow(I, [], 'Colormap', hot, 'initialMagnification', 1000);
%   bwMask = fhroi(f, gca);
%
% See also IMFREEHAND, CREATEMASK

% v0.6 called fhroi (Feb 2012) by Andrew Davis -- addavis@gmail.com -
% v0.7 (Feb 2018) modified and renamed as multiROI by Andrew Lin - andrew.lin@sheffield.ac.uk
% instead of using a menu, use keypresses to say add or subtract


% Check input and set up variables
if ~exist('axHandle', 'var')
    axHandle = gca; 
end
imHandle = imhandles(axHandle);
imHandle = imHandle(1);             % First image on the axes
hOVM = [];                          % no overlay mask yet
axis equal, axis tight;

% Create and activate a pointer manager for the figure:
iptPointerManager(figHandle, 'enable');
% Have the pointer change to a cross when the mouse enters an axes object:
iptSetPointerBehavior(axHandle, @(figHandle, currentPoint)set(figHandle, 'Pointer', 'cross'));

imageDim = [range(axHandle.YLim), range(axHandle.XLim)];
if ~isempty(varargin) && ~isempty(varargin{1})
    ROIs = varargin{1};
    if length(size(ROIs)) > 3
        size(ROIs)
        error('ROI input wrong # dimensions')
    end
    cumulative = logical(sum(ROIs,3));
    hOVM = alphamask(cumulative, [1 0 0], 0.6, axHandle);
else
    ROIs = [];
    cumulative = zeros(imageDim);
end

% User instructions and initial area
disp('1. Use zoom and pan tools if desired');
disp('2. Make sure no tools are selected')
disp('3. Left click and drag to add closed loop');
disp('4. Press a to add ROIs; press s to subtract ROIs; press t to finish');

set(figHandle, 'WindowKeyPressFcn', @keyPressCallback)
while true
    waitforbuttonpress
    clickedPt = axHandle.CurrentPoint; % click on a point
    center = [clickedPt(1,1),clickedPt(1,2)];
    circle = createCirclesMask(imageDim, center, circleSize);
%     size(cumulative)
%     size(circle)
    new = circle & ~cumulative; %what is new in the circle that didn't exist already
    cumulative = cumulative | new;
    delete(hOVM);
    hOVM = alphamask(cumulative, [1 0 0], 0.6, axHandle);
    uiwait(figHandle);
    lastCharacter = get(figHandle,'CurrentCharacter');
    switch lastCharacter
        case 's' % cancel (subtract)
            cumulative = cumulative & ~new; % remove the new pixels
            delete(hOVM);
            hOVM = alphamask(cumulative, [1 0 0], 0.6, axHandle);
        case 'a' % yes (add)
%             size(ROIs)
%             size(circle)
            ROIs = cat(3, ROIs, circle);
        case 't' % finished
            ROIs = cat(3, ROIs, circle);
            break;
    end
end

end

function keyPressCallback(source,eventdata)
% determine the key that was pressed
keyPressed = eventdata.Key;
uiresume(source);
switch keyPressed
    case 'a'
        disp('circle confirmed');
        
    case 's'
        disp('circle cancelled');
end
end