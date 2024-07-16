function [result, montageDims] = stackToMontage(stack, varargin)

% STACKTOMONTAGE: Turn a 3D stack into a montage
% Behaves the same as Matlab's builtin montage function in
% the mode montage(I) except it doesn't display the function
% 
% Usage:
% [montage, montageDims] = stackToMontage(stack, numRows, numCols)
%   stack is the input data and should be M x N x 1 x K for grayscale, and
%   M x N x 3 x K for colour
%   numRows is how many rows you want in the montage
%   numCols is how many columns you want
% [montage, montageDims] = stackToMontage(stack)
%   as above but the function will automatically determine the number of
%   rows and columns
% [montage, montageDims] = stackToMontage(stack, numRows)
%   numRows rows and ceil(numSlices/m) columns
%
% Outputs:
% montage is the montage
% montageDims is a size-3 array with the dimensions of how the montage was
% assembled
% montageDims(1) = numRows (how many tiles vertically)
% montageDims(2) = numCols (how many tiles horizontally)
% montageDims(3) = zDim
% 
% Usage notes:
% Save montageDims if you want to convert the montage back into a stack, or
% if you want to convert some manipulation of the montage (eg drawing an
% ROI) into a stack with the same dimensions as the original
% example usage:
% [montage, montageDims] = stackToMontage(stack);
% 
% If your stack is a 3D grayscale matrix, use the permute function to move
% the z-dimension to the 4th dimension:
% result = stackToMontage(permute(stack, [1 2 4 3]));
% 
% Matlab treats the first dimension as row indices, even though in images the
% first dimension is horizontal. To display the image in the correct
% orientation, use permute to switch dimension 1 and 2 (equivalent to the '
% transpose function which only works on 2D matrices)
% result = stackToMontage(permute(stack, [2 1 4 3]));

% Written by Andrew Lin


mDim = size(stack,1);
nDim = size(stack,2);
cDim = size(stack,3);
if (cDim~=1)&&(cDim~=3)
    error('stack must be M x N x 1 x K for grayscale or M x N x 3 x K for colour');
end
zDim = size(stack,4);

if isempty(varargin)
    % automatically calculate number of rows and columns
    % I started to try to do this more elegantly but gave up
    % refRows and refCols define how many rows and columns you want for up to 20 elements. After that, just use ceil(sqrt())
    refRows = [1,1,1,2,2,2,2,2,3,2,2,2,3,3,3,4,3,3,4,4];
    refCols = [1,2,3,2,3,3,4,4,3,5,6,6,5,5,5,4,6,6,5,5];
    if (zDim<=20)
        numRows = refRows(zDim);
        numCols = refCols(zDim);
    else
        numRows = ceil(sqrt(zDim));
        numCols = ceil(sqrt(zDim));
    end
elseif length(varargin)==1
    numRows = varargin{1};
    numCols = ceil(zDim/numRows);
elseif length(varargin)==2
    numRows = varargin{1};
    numCols = varargin{2};
else
    error('stackToMontage: wrong number of input arguments');
end


result = zeros(mDim*numRows, nDim*numCols, cDim);

for i=1:zDim
    mPos = floor((i-1)/numCols); %starts at 0
    nPos = mod((i-1),numCols); %starts at 0
    result((1+mPos*mDim):(mDim*(mPos+1)), (1+nPos*nDim):(nDim*(nPos+1)), :) = stack(:,:,:,i);
end

% get rid of the singleton dimension for grayscale images
result = squeeze(result);
montageDims = [numRows numCols zDim];

end