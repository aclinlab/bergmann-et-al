function [rgbimage] = overlayColorOnGrayscale(colorimage, grayimage, colorrange, grayrange)
% overlayColorOnGrayscale: overlay false-coloring (eg activity) on a grayscale image (eg baseline fluorescence)
% colorimage and grayimage are 3d matrices
% usage:
% rgbimage = overlayColorOnGrayscale(colorimage, grayimage, colorrange, grayrange)
% Revised Aug 2017 to handle 3d input. previous version only handled 2d
% matrices
% There is no reason to make this work for 4 inputs (ie adding different
% channels) because different channels will have different ranges. Just
% call this function multiple times for each channel.
% also revised to auto-contrast the images if you put scalar values for
% colorrange and grayrange
% resulting image is an rgb matrix (i.e. height x width x 3 x depth)
% where the background is the grayimage in grayscale
% and colorimage is overlaid on top
% colorrange and grayrange are e.g. [0 255] or [0 100] to define
% the min and max of colorimage and grayimage for the overlay image
% you can use this to enhance contrast or brighten the image, for example
% Written by Andrew Lin


if size(colorimage)~=size(grayimage)
    error('Two images are not equal size!');
end

numcolorpoints = 256;

maxColorProportion = 0.99;
maxGrayProportion = 0.99;

if isscalar(colorrange)
    sortColorImageArray = sort(colorimage(:));
    maxColor = sortColorImageArray(ceil(maxColorProportion*length(sortColorImageArray)));
    colorrange = [0 maxColor];
end
if isscalar(grayrange)
    sortGrayImageArray = sort(grayimage(:));
    maxGray = sortGrayImageArray(ceil(maxGrayProportion*length(sortGrayImageArray)));
    grayrange = [0 maxGray];
end
colorrange

% change 'hot' to the color map of your choice
color_map = hot(numcolorpoints);

height = size(colorimage,1);
width = size(colorimage,2);
depth = size(colorimage,3);

rgbimage = zeros(height, width, 3, depth);

colorimage(colorimage > colorrange(2)) = colorrange(2);
colorimage = colorimage - colorrange(1);
colorimage(colorimage < 0) = 0;
colorimage = colorimage .* (numcolorpoints / (colorrange(2) - colorrange(1)));

grayimage(grayimage > grayrange(2)) = grayrange(2);
grayimage = grayimage - grayrange(1);
grayimage = grayimage ./ (grayrange(2) - grayrange(1));

for i=1:height
    for j=1:width
        for k=1:depth
            if colorimage(i,j,k)==0
                rgbimage(i,j,1,k) = grayimage(i,j,k);
                rgbimage(i,j,2,k) = grayimage(i,j,k);
                rgbimage(i,j,3,k) = grayimage(i,j,k);
            else

                rgb = color_map(ceil(colorimage(i,j,k)),:);

                rgbimage(i,j,:,k) = rgb;
            end
        end
    end
end

maxvalue = max(rgbimage(:));
rgbimage = rgbimage./maxvalue;
rgbimage = squeeze(rgbimage); % if there was no 3rd dimension, get rid of it

end