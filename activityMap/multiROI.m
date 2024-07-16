function bwMask = multiROI(figHandle, axHandle)
% multiROI:  Interactively specify 2D freehand ROI
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

% User instructions and initial area
disp('1. Use zoom and pan tools if desired');
disp('2. Make sure no tools are selected')
disp('3. Left click and drag to add closed loop');
disp('4. Press a to add ROIs; press s to subtract ROIs; press z to undo the last drawing; press t to finish');

fhObj = drawfreehand(axHandle);          % choose initial area

% if you de-comment this, the program will wait and you can drag the
% outline around - when you are happy, double click to resume
% position = wait(fhObj);               % allow repositioning

set(figHandle, 'WindowKeyPressFcn', @keyPressCallback)
% (x,y)disp(get(fhObj,'CurrentCharacter'))
try
   bwMask = createMask(fhObj, imHandle);  % logical matrix of image size
catch ME
   error('bwMask ROI was not created properly');
end
delete(fhObj);                         % clean up
lastMask = bwMask;
% global RETURN_PRESSED
% RETURN_PRESSED = false;
% % Await user input to determine if the ROI needs tweaking
% while ~RETURN_PRESSED
% %     uiwait(figHandle);
while true
    delete(hOVM);                       % delete old overlay mask
    hOVM = alphamask(bwMask);           % overlay image with mask
    disp('Waiting for you to press a, s, z or t');
    uiwait(figHandle);
    fhObj = [];
    lastCharacter = get(figHandle,'CurrentCharacter');
    switch lastCharacter
        case {'z'}
            if bwMask==lastMask
                disp('Out of undos');
            else
                bwMask = lastMask;
            end
        case {'s','a'}
            fhObj = drawfreehand;
            try
                newMask = createMask(fhObj, imHandle);  % logical matrix of image size
            catch ME
                disp('bwMask ROI was not created properly');
                disp('Press a to add ROIs; press s to subtract ROIs; press z to undo the last drawing; press t to finish');
            end
            delete(fhObj);
            switch lastCharacter
                case 's'
                    lastMask = bwMask;
                    bwMask = bwMask & ~newMask;
                case 'a'
                    lastMask = bwMask;
                    bwMask = bwMask | newMask;
            end
%             fhSub = imfreehand;
%             bwSub = createMask(fhSub, imHandle);
%             bwMask = bwMask & ~bwSub;        % logical 'bwMask and not bwSub'
%             delete(fhSub);
%         case 'a'
%             fhAdd = imfreehand;
%             bwAdd = createMask(fhAdd, imHandle);
%             bwMask = bwMask | bwAdd;         % logical 'bwMask or bwAdd'
%             delete(fhAdd);
        case {'t'} % 't' for 'terminate
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
        disp('Now adding ROIs');
        
    case 's'
        disp('Now subtracting ROIs');
    case 'z'
        disp('You pressed undo');
end
end