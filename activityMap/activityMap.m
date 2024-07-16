classdef activityMap
    % activityMap: Class for calculating and storing dF/F activity for 5D movies
    % Calculate an activity map for a 5-dimensional ScanImage movie where time is
    % the last dimension (order of dimensions x y z c t)
    % 
    % Example usage:
    % map = activityMap(params); 
    % figure, imagesc(map.overlayMontage([0 3], [1 2]);
    % map.saveMap();
    % map.useOtherStdevThres(multiple); 
    % mask = activityMap.outlineObject(map.f0, channel);
    %
    % There is the option to calculate the activity map in a "control" stimulus period
    % (e.g. "empty bottle") and subtract this from the "real" stimulus
    % period
    % 
    % Key functions for the user (type 'help activityMap.xxxx' to see documentation for
    % function xxxx)
    % - activityMap (constructor)
    % - overlayMontage (view the activity map overlaid on baseline fluorescence)
    %
    
    % Written by Andrew Lin based on code inherited from Ricardo Marques.
    % Modified from what was used in Lin et al 2014 Nat Neurosci by adding the
    % ability to deal with 5-dimensional movies (including multiple z-slices or
    % channels)
    properties  
        % key data that's calculated
        
        f0 % average signal during prePeriod
        f0unfiltered % average signal during prePeriod without xy smoothing
        f0stdev % stdev (over time) of signal during prePeriod
        dff % dF/F during stimPeriod (minus dF/F during ctrlPeriod, if applicable)
        dffNoThres % dF/F with no threshold applied
        responsivePixelsMeanTimecourse % mean dF/F over time of all pixels with dff > 0
        
        % parameters of the calculation
        
        correlations % correlations of each frame with f0.
        mismatches % the frames being discarded
        prePeriod % the frame numbers for the pre-stimulus period (defined by params.preLimits)
        stimPeriod % the frame numbers for the stimulus period (defined by params.stimLimits)
        ctrlPeriod % the frame numbers for the control stimulus period (defined by params.ctrlLimits)
        numStimPeriods % number of stim periods
        
        % movie parameters  
        
        dimensions % dimensions of the movie, usually x, y, z, c, t
        frameRate % frames per s (Hz), taken from the metadata of the ScanImage file
        nRows % length on x dimension (Matlab treats this as rows so you need to flip rows and columns when displaying the image
        nColumns % length on y dimension
        nSlices % z dimension
        nChannels % number of channels
        nFrames % number of frames in the whole movie
        rowCal % calibration (microns per pixel) on x dimension (Matlab treats this as rows so you need to flip rows and columns when displaying the image
        columnCal % calibration (microns per pixel) on y dimension
        sliceCal % calibration (microns per pixel) on z dimension
        

        % parameter object - of class activityMapParams
        params %@activityMapParams
        % which has properties:
        %         % these need to be set before constructing an activityMap
        %         preLimits@double
        %         stimLimits@double
        %         ctrlLimits@double % optional
        %
        %         % these have default values but you might want to set them
        %         corrThres@double = 0
        %         stdevThresMult@double = 2
        %         bkgnd@double = 0
        %         capForCorrCalc@double = 10
        %         filter = fspecial('gaussian',5,2)
        %         userDrawBkgnd@logical = false
        %         userDrawBkgndChannel = 1
        %
        %         % these can be set by the activityMap constructor function
        %         fileName@char
        %         pathName@char
        %         userBkgndROI
        %    
        %         % immutable
        %         softwareVersion
        %  
        % optional properties for some types of analysis
        
        skel % skeleton object for skeleton analysis
        voronoiTimeSeries % dF/F of each Voronoi space defined by spaced nodes in the user-defined skeleton
        voronoiMeanResp % mean resp of each Voronoi space during stimPeriod
        voronoiMaxResp % max resp of each Voronoi space during stimPeriod
        offsets = [0 0]
        ymin= 1
        ymax
        xmin= 1
        xmax
        
        ROIs
        ROITimeSeries
        ROIStimPeriod
        intraROIcorr
        
    end
    
    properties (Dependent)
        pixelCal % [rowCal columnCal sliceCal]
    end
    
    methods

        % to do: need a way to exclude the slices and pixels that are empty
        % due to the xy- or z-registration. Problem, this is not included
        % in the metadata! Maybe the registration algorithms should track
        % this and delete the relevant pixels / slices?
        % also to do: bkgnd should be a vector that's numChannels long
        function obj = activityMap(params)
            % Takes an activityMapParams object as input
            % activityMap:
            % Constructor for the class activityMap.
            % Takes an activityMapParams object as input
            % if there are no input arguments, creates a completely empty
            % activityMap object by default
            
            % default: return empty object
            if ~nargin
                return;
            end
            % initialise parameters
            obj.params = params;
            
            % error check the input parameters
            if isempty(obj.params.preLimits)||isempty(obj.params.stimLimits)
                error('F0 or stimulus time limits undefined');
            end
            
            % if filepath is undefined in input parameters, ask user to
            % select file
            if isempty(params.fileName) || isempty(params.pathName)
                [obj.params.fileName, obj.params.pathName] = uigetfile( ...
                    '*.tif', 'Pick a .tif file');
            end
            
            if isempty(params.preLimits) || isempty(params.stimLimits)
                error('F0 or stimulus period undefined in input parameters');
            end
                
            % read the movie
            fprintf(strcat('Loading file ', obj.params.pathName, obj.params.fileName, '...\n'));
            [movie, SI] = readScanImageTiffLinLab(strcat(obj.params.pathName,obj.params.fileName));
            % readScanImageTiffLinLab returns a 5D matrix even if some of
            % these dimensions are singletons
            % It is crucial to keep the dimension order rigid and not
            % squeeze out singleton dimensions. What if you have only 1
            % z-slice and 2 channels, vs 2 z-slices and 1 channel. If you
            % squeeze out the singleton dimension, how will you know which
            % is the case?
            movie = single(movie);
            obj.frameRate = SI.hRoiManager.scanVolumeRate;
            
            
            % Get dimensions of the movie
            obj.dimensions = size(movie);
            numDims = length(obj.dimensions);
            % time is the last dimension
            timeDimension = numDims;
            meanOverTime = mean(movie,timeDimension);
            obj.nRows = size(meanOverTime,1);
            obj.nColumns = size(meanOverTime,2); % length on y dimension
            obj.nSlices = size(meanOverTime,3);% z dimension
            obj.nChannels = size(meanOverTime,4);% number of channels
            obj.ymax = obj.nRows;
            obj.xmax = obj.nColumns;
            
            % get x y z calibration
            if obj.nSlices>1
                obj.sliceCal = SI.hStackManager.stackZStepSize;
            end
            % x and y calibration based on measurements of 170929
            % at 8x zoom, 256 x 256 pixels, 0.32 microns per pixel
            % this is the first dimension which is x or fast scan (but
            % Matlab treats it as rows)
            obj.rowCal = 0.32 * SI.hRoiManager.scanAngleMultiplierFast / (SI.hRoiManager.scanZoomFactor/8) / (SI.hRoiManager.pixelsPerLine/256);
            % this is the second dimension which is y or slow scan (but
            % Matlab treats it as columns)
            obj.columnCal = 0.32 * SI.hRoiManager.scanAngleMultiplierSlow / (SI.hRoiManager.scanZoomFactor/8) / (SI.hRoiManager.linesPerFrame/256);

            obj.nFrames = obj.dimensions(timeDimension);
            
            % get the frame numbers corresponding to the given limits for
            % pre and stim
            obj.prePeriod = obj.getFramesFromLimits(obj.params.preLimits);
            if iscell(obj.params.stimLimits)
                obj.numStimPeriods = length(obj.params.stimLimits);
                obj.stimPeriod = cell(1,obj.numStimPeriods);
                for i=1:obj.numStimPeriods
                    obj.stimPeriod{i} = obj.getFramesFromLimits(obj.params.stimLimits{i});
                end
            else
                obj.numStimPeriods = 1;
                obj.stimPeriod = obj.getFramesFromLimits(obj.params.stimLimits);
            end
                      
            % ********************************************************
            % Begin processing
            % ********************************************************
            
            fprintf('Filtering...\n');
            tic
            % Filter the movie throughout, then take fpre as an interval of this
            % NB h1 is a 2D filter so imfilter only operates on the first 2 dimensions
            % (should be x and y) - doesn't smooth across z, channels, or time
            filtered = imfilter(movie,obj.params.filter,'replicate');
            unfiltered = movie;
            
            % smooth data in time dimension (5th dimension) over 5 frames
            if isfield(obj.params,'timeSmoothingWindow') && (obj.params.timeSmoothingWindow > 1)
                filtered = smoothdata(filtered, 5, 'movmean', obj.params.timeSmoothingWindow);
                unfiltered = smoothdata(unfiltered, 5, 'movmean', obj.params.timeSmoothingWindow);
            end
            
            toc
            
            
            if ~isempty(obj.params.userBkgndROI)
                for i=1:obj.nChannels
                    thisChannel = meanOverTime(:,:,:,i); % if meanOverTime is only 3D but i is 1, this is fine
                    bkgndPixels = thisChannel(obj.params.userBkgndROI==1);
                    obj.params.bkgnd(i) = mean(bkgndPixels(:));
                end
            elseif (obj.params.userDrawBkgnd)
                % if there are multiple channels, use the channel defined
                % in the parameters
                if obj.nChannels>1
                    forMontage = meanOverTime(:,:,:,obj.params.userDrawBkgndChannel);
                else
                    forMontage = squeeze(meanOverTime);
                end 
                 
                [montageForBkgnd, montageDims] = stackToMontage(permute(forMontage, [2 1 4 3]));
%                 montageDims
                fig = figure('Name','Draw background ROI');
                imagesc(montageForBkgnd);
                axis equal, axis tight;
                
                % ask user to select area
                maskMontage = multiROI(fig, gca);
                close(fig);
%                 size(maskMontage)
                
                % ROI will be 3D
                obj.params.userBkgndROI = permute(montageToStack(maskMontage, ...
                    montageDims), [2 1 3]);
                size(obj.params.userBkgndROI)
                figure,imagesc(stackToMontage(permute(obj.params.userBkgndROI,[2 1 4 3])));
                for i=1:obj.nChannels
                    thisChannel = meanOverTime(:,:,:,i); % if meanOverTime is only 3D but i is 1, this is fine
                    bkgndPixels = thisChannel(obj.params.userBkgndROI==1);
                    obj.params.bkgnd(i) = mean(bkgndPixels(:));
                end
            end
            
            fprintf('Subtracting background...\n');
            tic
            % Subtract background AFTER filtering - important for images with high
            % background, as this high background is usually noisy, so you end up
            % with too many negative pixels
            if isscalar(obj.params.bkgnd)
                filtered = filtered - obj.params.bkgnd;
                unfiltered = unfiltered - obj.params.bkgnd;
            else
                % a vectorized implementation of the for loop below
                % filtered = filtered - repmat(permute(obj.params.bkgnd(:),[2 3 4 1]),obj.nRows,obj.nColumns,obj.nSlices,1);
                for i=1:obj.nChannels
                    filtered(:,:,:,i) = filtered(:,:,:,i) - obj.params.bkgnd(i);
                    unfiltered(:,:,:,i) = unfiltered(:,:,:,i) - obj.params.bkgnd(i);
                end
            end
            toc
            % in this section we remove slices that deviate excessively
            % from the mean F0 image - that is, slices for which the
            % correlation with the mean F0 image falls below
            % params.corrThres
            
            fprintf('Calculating F0, inter-frame correlations...\n');
            keepCheckingMismatches = 1;
            while keepCheckingMismatches
                % the {':'} notation is to handle inputs with any number of
                % dimensions
                index = repmat({':'},1,timeDimension);
                index{timeDimension} = obj.prePeriod;
                fpre = filtered(index{:});
                fpreunfiltered = unfiltered(index{:});
                
                % Take the mean of the pre interval along the time dimension
                % This results in a  matrix as a "snapshot" of the movie during pre
                obj.f0 = mean(fpre,timeDimension);
                obj.f0unfiltered = mean(fpreunfiltered,timeDimension);
                % Don't zero out the raw filtered data, only zero out the f0
                obj.f0(obj.f0 < 0) = 0;
                obj.f0unfiltered(obj.f0unfiltered < 0) = 0;
                
                obj.correlations = obj.corrOverTime(filtered, obj.f0, obj.params.capForCorrCalc);
                obj.mismatches = find(obj.correlations < obj.params.corrThres);
                
                newprePeriod = obj.prePeriod;
                for i=1:length(obj.mismatches)
                    newprePeriod(newprePeriod==obj.mismatches(i)) = [];
                end
                
                if (size(newprePeriod)==size(obj.prePeriod))
                    keepCheckingMismatches = 0;
                end
                
                obj.prePeriod = newprePeriod;
            end
            
            if iscell(obj.stimPeriod)
                removeFromStim = intersect(obj.mismatches,cell2mat(obj.stimPeriod));
                if ~isempty(removeFromStim)
                    for i=1:obj.numStimPeriods
                        obj.stimPeriod{i} = setdiff(obj.stimPeriod{i},removeFromStim');
                    end
                end
            else
                removeFromStim = intersect(obj.mismatches, obj.stimPeriod);
                if ~isempty(removeFromStim)
                    obj.stimPeriod = diff(obj.stimPeriod,removeFromStim);
                end
            end
            if ~isempty(removeFromStim)
                fprintf('Warning, removed these frames from stim period:\n');
                removeFromStim; %intersect(obj.mismatches,obj.stimPeriod)
            end
%             for i=1:length(obj.mismatches)
%                 obj.stimPeriod(obj.stimPeriod==obj.mismatches(i)) = [];
%             end
            toc
            fprintf('Calculating stdev...\n');
            tic
            % Take the standard deviation of the prestimulus interval along the time
            % dimension
            obj.f0stdev = std(fpre,1,timeDimension);
            toc
            fprintf('Calculating DF...\n');
            tic
            % Calculate DF
            indexStimPeriod = repmat({':'},1,numDims);
            df = zeros(obj.nRows, obj.nColumns, obj.nSlices, obj.nChannels, obj.numStimPeriods);
            if iscell(obj.stimPeriod)
                for i=1:obj.numStimPeriods
                    indexStimPeriod{numDims} = obj.stimPeriod{i};
                    df(:,:,:,:,i) = mean(filtered(indexStimPeriod{:}),timeDimension) - obj.f0;
                end
            else
                indexStimPeriod{numDims} = obj.stimPeriod;
                df = mean(filtered(indexStimPeriod{:}),timeDimension) - obj.f0;
            end

            % if there is a defined 'control stimulus' period, get the dF
            % for that time
            if ~isempty(params.ctrlLimits)
                obj.ctrlPeriod = obj.getFramesFromLimits(obj.params.ctrlLimits);
                for i=1:length(obj.mismatches)
                    obj.ctrlPeriod(obj.ctrlPeriod==obj.mismatches(i)) = [];
                end
                indexCtrlPeriod = repmat({':'},1,numDims);
                indexCtrlPeriod{numDims} = obj.ctrlPeriod;
                ctrldf = mean(filtered(indexCtrlPeriod{:}),timeDimension) - obj.f0;
                % zero out negative pixels
                ctrldf(ctrldf<0) = 0;
            else
                ctrldf = 0;
            end
            
            df = df - repmat(ctrldf,1,1,1,1,obj.numStimPeriods);
            % divide by F0 to get df/f
            obj.dff = df./repmat(obj.f0,1,1,1,1,obj.numStimPeriods);
            obj.dffNoThres = obj.dff;
            
            % Zero any pixels where the signal is less than stdevThresMult*stdev in the pre period
            stdThres = obj.params.stdevThresMult*obj.f0stdev./obj.f0;
            obj.dff(obj.dff<repmat(stdThres,1,1,1,1,obj.numStimPeriods)) = 0;
            obj.dff(isnan(obj.dff)) = 0;
            obj.dff(isinf(obj.dff)) = 0;
            toc
            
%             responsivePixels = find(obj.dff > 0);
%             [xs,ys,zs,cs] = ind2sub(size(obj.dff), responsivePixels);
%             
%             nResponsivePixels = length(responsivePixels);
%             timeCourses = zeros(nResponsivePixels,obj.nFrames,obj.nChannels);
%             for i=1:nResponsivePixels
%                 % calculate dff as (F-F0)/F0 over whole time course
%                 for j=1:obj.nChannels
%                     timeCourses(i,:,j) = (filtered(xs(i),ys(i),zs(i),j,:) - obj.f0(xs(i),ys(i),zs(i),j))/ ...
%                         obj.f0(xs(i),ys(i),zs(i),j);
%                 end
%             end
%             obj.responsivePixelsMeanTimecourse = mean(timeCourses,1);
        end
        
        function obj = subPlanes(obj, z)
            % return a version of the activity map at only a certain planes
            % given by z
            % e.g.
            % mapSinglePlane = map.subPlanes(2); % return only z-plane 2
            % mapSinglePlane = map.subPlanes([2 4 5]); % planes 2, 4, 5
            obj.f0 = obj.f0(:,:,z);
            obj.f0stdev = obj.f0stdev(:,:,z);
            obj.dff = obj.dff(:,:,z);
            obj.dffNoThres = obj.dffNoThres(:,:,z);
            obj.dimensions(3) = length(z);
            obj.nSlices = length(z);
        end
        
        function result = getFramesFromLimits(obj,limits)
            % return frame numbers between limits (in seconds)
            % getFramesFromLimits:
            % limits is an nx2 matrix
            % each row is the beginning and end of a time period, in seconds
            % this allows your time period to span multiple discontinuous
            % periods
            % eg if the input is [5 10; 12 17; 30 34]
            % then this function will return frames that occur between 5-10 s,
            % 12-17 s, and 30-34 s
            limitsDim = size(limits);
            result = [];
            if nnz(isnan(limits))>0
                return;
            end
            if max(limits(:)) > obj.nFrames*obj.frameRate
                error('period limits exceed total length of movie');
            end
            for i=1:limitsDim(1)
                result = [result ...
                    round(limits(i,1)*obj.frameRate):round(limits(i,2)*obj.frameRate)];
            end
            % double-check: remove any frames that occur before the 
            % beginning or end of the movie
            result(result<1) = [];
            result(result>obj.nFrames) = [];
        end
        
        function result = useOtherStdevThres(obj, otherStdevThresMult)
            % returns dF/F with a different standard deviation threshold
            % Usage:
            % newDFF = map.useOtherStdevThres(newMultiple);
            %  default is that during the stimulus period, pixels must cross 2x the
            %  standard deviation during the pre-stimulus period. This could have
            %  been set differently in the constructor parameters. Using this
            %  method, change the threshold to 'multiple' instead of 2x

            result = obj.dffNoThres;
            result(result<(otherStdevThresMult*obj.f0stdev./obj.f0)) = 0;
        end
        

        function [result, zMontageDims, channelMontageDims] = overlayMontage(obj, varargin)
            % Returns montage of activity map with DF/F overlaid on baseline fluorescence
            % Usage:
            % obj.overlayMontage(): returns montage with all channels and
            % automatically determined display ranges
            %
            % obj.overlayMontage(displayRange): returns montage with all
            % channels where the df/f is set to displayRange, which is a
            % numChannels x 2 matrix, so that each channel has its own
            % display range. Each channel's df/f is overlaid on its own
            % baseline fluorescence
            %
            % obj.overlayMontage(displayRange, channels): as above except
            % displayRange is an n x 2 matrix and 
            % channels is an n x 2 matrix for each row is a pair of channel
            % indices, the first of which is the channel for which you want
            % to display the df/f, and the second is the channel for which
            % you want to display the f0
            % 
            % To preserve the dimensions of the automatically-generated
            % montages, call the function like this:
            % [result, zMontageDims, channelMontageDims] = ...
            %       obj.overlayMontage(); % or with optional input parameters
            %           
            % 
            % Example:
            % result = obj.overlayMontage([0 8; 0 2], [1 2; 2 2]);
            % would give you:
            % df/f of channel 1 overlaid on f0 of channel2, with a display
            % range of [0 8] for the df/f
            % and df/f of channel 2 overlaid on f0 of channel2, with a
            % display range of [0 2] for the df/f
            %
            % Each row in displayRange corresponds to the matching row in
            % channels
            %
            % To display, run e.g.:
            % figure, imagesc(obj.overlayMontage([0 8; 0 2], [1 2; 2 2]));
            numDefaultChannels = size(obj.dff,4);
            mapNum = 1;
            switch length(varargin)
                case 0
                    % each channel will have a displayRange of 0 because
                    % overlayColorOnGrayscale takes scalar 0 to mean set
                    % the display range automatically
                    displayRange = zeros(numDefaultChannels,1);
                    % make channels something like:
                    % [1 1
                    %  2 2];
                    channels = repmat((1:numDefaultChannels)', 1, 2);
                case 1
                    displayRange = varargin{1};
                    channels = repmat((1:numDefaultChannels)', 1, 2);
                case 2
                    displayRange = varargin{1};
                    channels = varargin{2};
                case 4
                    displayRange = varargin{1};
                    if ~isempty(varargin{2})
                        channels = varargin{2};
                    else
                        channels = repmat((1:numDefaultChannels)', 1, 2);
                    end
                    mapNum = varargin{3};
                    montageRows = varargin{4}(1);
                    montageCols = varargin{4}(2);
                    
            end
            numChannels = size(channels,1);
            for i=1:numChannels
                overlayImage = overlayColorOnGrayscale(permute(obj.dff(:,:,:,channels(i,1),mapNum), [2 1 3 4]), ...
                    permute(obj.f0(:,:,:,channels(i,2)), [2 1 3 4]), displayRange(i,:), 0);
                if length(varargin)==4
                    [overlayStack(:,:,:,i), zMontageDims] = stackToMontage(overlayImage,montageRows,montageCols);
                else
                    [overlayStack(:,:,:,i), zMontageDims] = stackToMontage(overlayImage);
                end
            end
            [result, channelMontageDims] = stackToMontage(overlayStack, numChannels);
        end
        
        function [] = displayMultiMaps(obj, range)
            figure
            numMaps = size(obj.dff,5);
            for i=1:numMaps
                subplot(numMaps,1,i);
                imagesc(obj.overlayMontage(range,[],i,[1 size(obj.dff,3)]));
                axis equal, axis tight;
                set(gca,'xtick',[]);
                set(gca,'ytick',[]);

            end
        end
        
        
        function saveMap(obj, varargin)
            % save this map as the name of the original file, minus the .tif extension, 
            % plus '-activityMap.mat' plus optional suffix eg
            % '-activityMap-2.mat')
            % Usage:
            % map.saveMap(); - use default suffix, '-activityMap.mat'
            % map.saveMap(suffix); - add the suffix between 'activityMap'
            % and '.mat'
            % eg map for 171023f1_00002.tif saved as
            % 171023f1_00002-activityMap.mat
            if nargin==2
                suffix = varargin{1};
            else
                suffix = '';
            end
            extensionIndex = strfind(obj.params.fileName,'.');
            extensionIndex = extensionIndex(end); % get the index of the *last* dot
            save(strcat(obj.params.fileName(1:(extensionIndex-1)), ...
                '-activitymap', suffix, '.mat'),'obj');
        end
        
        function saveTiffStack(obj, colorRange, grayRange)
            % saveTiffStack
            % This function will save a .tif stack of the activity map
            % (dff) in "hot" false coloring overlaid on the baseline
            % fluorescence (f0) in grayscale
            % Usage:
            % map.saveTiffStack(colorRange, grayRange)
            % colorRange and grayRange are 1 x 2 arrays that define the minimum
            % and maximum of of dff (colorRange) or f0 (grayRange).
            % This is useful for making figures so you can force all the
            % panels in your figure to follow the same false coloring
            % scale.
            % Note: grayscale max values will depend on if the movie was
            % imported as 8-bit or 16-bit (16-bit pixel values are bigger)
            overlayImage = overlayColorOnGrayscale(permute(obj.dff, [2 1 3 4]), ...
                permute(obj.f0, [2 1 3 4]), colorRange, grayRange);
            filePrefix = obj.params.fileName(1:end-4);
            imwrite(overlayImage(:,:,:,1), strcat(filePrefix,'-overlayMap.tif'), 'tif', ...
                'Compression', 'lzw');
            for i=2:size(overlayImage,4)
                imwrite(overlayImage(:,:,:,i), strcat(filePrefix,'-overlayMap.tif'), 'tif', ...
                    'Compression', 'lzw','WriteMode', 'append');
                
            end
        end
        
        function plotCorrelations(obj)
            % plots the correlations for each frame vs F0
            % 2 x-axes: main one is time (s) and secondary is frame number
            % Use this to see if there are any sudden movements that you
            % might be able to exclude by setting the corrThres parameter
            % (ie throw away frames that have less than params.corrThres
            % correlation with the F0)
            frameSeries = 1:obj.nFrames;
            timeSeries = frameSeries/obj.frameRate;
            figure;
            
            [frameAxis, timeAxis] = obj.dualXaxes();
            plot(timeSeries,obj.correlations);
            ylabel(timeAxis,'Correlation of each frame to F0');
           
%             % extra axis for frame number
%             frameAxis = axes('Position', [.1 .1 .8 1e-12]);
%             set(frameAxis,'Units','normalized');
%             xlabel(frameAxis,'Frame number');
%             xlim([0 obj.nFrames]);
%             
%             % time axis
%             timeAxis = axes('Position', [.1 .2 .8 .7]);
%             set(timeAxis, 'Units', 'normalized');
%             plot(timeSeries,obj.correlations);
%             xlabel(timeAxis,'Time (s)');
%             ylabel(timeAxis,'Correlation of each frame to F0');
            
        end
        
        function [frameAxis, timeAxis] = dualXaxes(obj)
            % Create dual time and frame axes for plotting time series data
            
            % extra axis for frame number
            frameAxis = axes('Position', [.1 .1 .8 1e-12]);
            set(frameAxis,'Units','normalized');
            %xlabel(frameAxis,'Frame number');
            %xlim(frameAxis, [0 obj.nFrames]);
            
            % time axis
            timeAxis = axes('Position', [.1 .2 .8 .7], 'FontSize',17);
            set(timeAxis, 'Units', 'normalized');
            xlabel(timeAxis,'Time (s)');
            %xlim(timeAxis, [0 obj.nFrames/obj.frameRate]);
        end
        

        
        function [] = displayOutline(obj, mask, channel)
            % display a mask on montage view
            % Usage:
            % map.displayOutline(mask, channel)
            %  mask: logical matrix (0s and 1s) - often stored in
            %  map.skel.mask
            %  channel: which channel's F0 do you want to display?
            
            channelForDisplay = obj.f0(:,:,:,channel);
            if size(channelForDisplay)~=size(mask)
                error('displayOutline: mask and f0 must be the same size');
            end
            [montageForDisplay, montageDims] = stackToMontage(permute(channelForDisplay,[2 1 4 3]));
            figure('Name', [obj.params.fileName, 'F0 channel ', num2str(channel)]);
            imagesc(montageForDisplay);
            axis equal, axis tight;
            maskMontage = stackToMontage(permute(mask,[2 1 4 3]), montageDims(1), montageDims(2));
            alphamask(maskMontage);
            title(['Mask for ', obj.params.fileName, ', F0 channel ', num2str(channel)]);
        end
        
        function obj = makeSkeleton(obj, channel, spacing, varargin)
            % Prompt user to define a skeleton (instructions will appear on
            % the console)
            % Usage:
            % obj = obj.makeSkeleton(channel, spacing)
            %  channel: which channel you want to use to draw the mask and
            %           starting point
            %  spacing: program will use the user-defined skeleton to
            %           create evenly spaced nodes. What spacing do you
            %           want?
            % Alternative usage:
            % obj = obj.makeSkeleton(channel, spacing, startingPointChannel)
            %  channel, spacing as above
            %  startingPointChannel - which channel do you want to use to
            %  define the starting point? - if different than the channel
            %  for making the mask
            % Example:
            % obj = obj.makeSkeleton(2, 20);
            % -->use channel 2, spacing = 20 microns
            % obj = obj.makeSkeleton(1, 20, 2);
            % -->use channel 1 to define the mask, spacing = 20 microns,
            % use channel 2 to define the starting point
            if length(varargin)==1
                startingPointChannel = varargin{1};
            elseif isempty(varargin)
                startingPointChannel = channel;
            end
            disp('Step 1: Draw a mask');
            mask = activityMap.outlineObject(obj.f0, channel);
            disp('Step 2: Click points to define the skeleton');
            obj.skel = skeleton(mask, obj.pixelCal);
            disp('Step 3: Define the stimulation site');
            obj.skel = obj.skel.defineStartingPoint(obj.f0(:,:,:,startingPointChannel));
            obj.skel = obj.skel.labelBranches();
            obj.skel = obj.skel.createSpacedNodes(spacing);
            obj.skel = obj.skel.createVoronoiMask();
            obj.skel.drawSkeleton();
        end
        
        function obj = skeletonAnalysis(obj, channel, varargin)
            % Analyze dF/F of evenly spaced divisions on the skeleton
            % skeletonAnalysis will calculate both mean and
            % max and store them in the map object, but it will call
            % displaySkeletonAnalysis with only the default, ie display
            % mean, not max. Call displaySkeletonAnalysis again if you want
            % to show max
            % Usage:
            % obj = obj.skeletonAnalysis(channel)
            %  channel: which channel you want to use for the dF/F analysis
            % obj = obj.skeletonAnalysis(channel, stimLimits)
            %  channel: which channel you want to use for the dF/F analysis
            %  optional second parameter: stimLimits, a 1x2 array if you
            %  want to use a different stimulus period than what's stored
            %  in the activityMap object
            %  ex:
            %  obj = obj.skeletonAnalysis(1, [19 20])
            %
            % We have to read the movie in again because it wasn't stored.
            % maybe a bit wasteful of time? but perhaps more wasteful of
            % memory to store the entire movie.
            
            if length(varargin)==1
                framesToUse = obj.getFramesFromLimits(varargin{1});
            else
                framesToUse = obj.stimPeriod;
            end
            
            [movie, ~] = readScanImageTiffLinLab(strcat(obj.params.pathName,obj.params.fileName));
            movie = single(movie);
            dffChannel = squeeze(movie(:,:,:,channel,:));
            dffChannel = dffChannel-obj.params.bkgnd(channel);
            % apply offsets in case the activityMap object has been
            % registered
            dffChannel = dffChannel((obj.ymin - obj.offsets(1)):(obj.ymax - obj.offsets(1)), ...
                (obj.xmin - obj.offsets(2)):(obj.xmax - obj.offsets(2)),:,:);
            % calculate dF/F for each division defined by the spaced nodes
            % in the skeleton
            nDivisions = size(obj.skel.voronoiMask,4);
            obj.voronoiTimeSeries = zeros(nDivisions,obj.nFrames);
            if size(obj.skel.voronoiMask(:,:,:,1))~=size(dffChannel(:,:,:,1))
                error('activityMap:skeletonAnalysis: size of skeleton mask does not match size of movie');
            end
            for t=1:obj.nFrames
                thisFrame = dffChannel(:,:,:,t);
                for x=1:nDivisions
                    obj.voronoiTimeSeries(x,t) = mean(thisFrame(obj.skel.voronoiMask(:,:,:,x)==1));
                end
            end
            
            voronoiF0 = mean(obj.voronoiTimeSeries(:,obj.prePeriod),2);
            
            for t = 1:obj.nFrames
                obj.voronoiTimeSeries(:,t) = (obj.voronoiTimeSeries(:,t) - voronoiF0)./voronoiF0;
            end
            % uncomment this if you want to smooth the data
            %obj.voronoiTimeSeries = smoothdata(obj.voronoiTimeSeries,2);
           
            % calculate mean dF/F during stimulus period
            obj.voronoiMeanResp = mean(obj.voronoiTimeSeries(:,framesToUse),2);
            % calculate max dF/F during stimulus period
            obj.voronoiMaxResp = max(obj.voronoiTimeSeries(:,framesToUse),[],2);
            
            %obj.displaySkeletonAnalysis();
        end
        
        function displaySkeletonAnalysis(obj, varargin)
            % display the skeleton analysis 
            % (use this if analysis was already calculated)
            % Usage:
            % obj.displaySkeletonAnalysis(); % will display the mean
            % response
            % obj.displaySkeletonAnalysis('max'); % will display the peak
            % response

            useMax = 0;
            if length(varargin)==1
                if strcmp(varargin{1},'max')
                    useMax = 1;
                end
            end
            if useMax
                dataToDisplay = obj.voronoiMaxResp;
            else
                dataToDisplay = obj.voronoiMeanResp;
            end

            nDivisions = size(obj.skel.voronoiMask,4);

            % plot individual dF/F traces against time
            figure
            frames = 1:obj.nFrames;
            timePoints = frames/obj.frameRate;
            [~, timeAxis] = obj.dualXaxes();
            ylabel(timeAxis,'dF/F');
            
            hold on;
            cmap = colormap;
            maxDist = max([obj.skel.spacedNodes.distFromStart]);

            for x=1:nDivisions
                color = cmap(round(obj.skel.spacedNodes(x).distFromStart*63/maxDist+1),:);
                plot(timePoints, obj.voronoiTimeSeries(x,:),'Color',color);
            end
            title(obj.params.fileName);
           
            % plot mean dF/F during stimPeriod against distance from start
            % site
            obj.graphDataVsDistance(dataToDisplay);
            % label y-axis of the plot produced by graphDataVsDistance
            if useMax
                ylabel('Maximum dF/F during user-defined period');
            else
                ylabel('Mean dF/F during user-defined period');
            end
            
            % draw skeleton color-coding the Voronoi spaces according to
            % meanResp
            obj.skel.drawSkeleton(dataToDisplay);
            
            % display a table on the console
            % every row will be a single Voronoi division
            % column 1: distance from start
            % column 2: mean dF/F
            % column 3: branch number
            [[obj.skel.spacedNodes.distFromStart]' dataToDisplay [obj.skel.spacedNodes.branchNumber]']
        end
        
        function [] = graphDataVsDistance(obj, data)
            % plot arbitrary data against distance from start site
            % usage:
            % obj.graphDataVsDistance(data)
            if length(data)~=length([obj.skel.spacedNodes.distFromStart])
                error('data must have same length as number of spaced nodes');
            end
            figure;
            title(obj.params.fileName);
            hold on;
            scatter([obj.skel.spacedNodes.distFromStart], data);
            % draw connecting lines according to the skeleton
            for i = 1:size(obj.skel.spacedLinks,1)
                x1 = obj.skel.spacedNodes(obj.skel.spacedLinks(i,1)).distFromStart;
                x2 = obj.skel.spacedNodes(obj.skel.spacedLinks(i,2)).distFromStart;
                y1 = data(obj.skel.spacedLinks(i,1));
                y2 = data(obj.skel.spacedLinks(i,2));
                line([x1 x2], [y1 y2]);
            end
            xlabel('Distance from start site (microns)');
        end
        
        function result = get.pixelCal(obj)
            % returns [obj.rowCal, obj.columnCal, obj.sliceCal]
            result = [obj.rowCal, obj.columnCal, obj.sliceCal];
        end
        
        function result = compareSkeletonAnalysis(obj, obj2, maxOrMean, ratioOrDiff)
            % compareSkeletonAnalysis
            % Usage:
            % obj.compareSkeletonAnalysis(obj2, maxOrMean, ratioOrDiff)
            % compares the results of the skeleton analysis for obj and
            % obj2
            % maxOrMean: 'max' if you want to use the peak dF/F response
            %            'mean' if you want to use the mean
            % ratioOrDiff: 'ratio' if you want to use the ratio obj2/obj
            %              'diff' if you want to use obj2 - obj
            if ~isequal(obj.skel, obj2.skel)
                error('The two activityMap objects must have the same skeleton');
            end
            switch maxOrMean
                case 'max'
                    data1 = obj.voronoiMaxResp;
                    data2 = obj2.voronoiMaxResp;
                otherwise
                    data1 = obj.voronoiMeanResp;
                    data2 = obj2.voronoiMeanResp;
            end
            switch ratioOrDiff
                case 'ratio'
                    result = data2./data1;
                otherwise
                    result = data2 - data1;
            end
            length(result)
            obj.graphDataVsDistance(result);
            ylabel({[ratioOrDiff ', ' maxOrMean]; obj.params.fileName; ' vs. '; obj2.params.fileName});
            obj.skel.drawSkeleton(result);    
        end
        
        function obj = offsetMap(obj,offsets, xmin, xmax, ymin, ymax)
            %Registers the activityMap object to a set of offsets
            %Example: map = map.offsetMap(offsets, xmin, xmax, ymin, ymax)
            obj.f0 = obj.f0((ymin - offsets(1)):(ymax - offsets(1)), ...
                (xmin - offsets(2)):(xmax - offsets(2)),:,:);
            obj.f0stdev = obj.f0stdev((ymin - offsets(1)):(ymax - offsets(1)), ...
                (xmin - offsets(2)):(xmax - offsets(2)),:,:);
            obj.dff = obj.dff((ymin - offsets(1)):(ymax - offsets(1)), ...
                (xmin - offsets(2)):(xmax - offsets(2)),:,:);
            obj.dffNoThres = obj.dffNoThres((ymin - offsets(1)):(ymax - offsets(1)), ...
                (xmin - offsets(2)):(xmax - offsets(2)),:,:);
            obj.dimensions(1:2)=[ymax-ymin+1 xmax-xmin+1];
            obj.nRows = obj.dimensions(1);
            obj.nColumns = obj.dimensions(2);
            obj.offsets = offsets;
            obj.xmin = xmin;
            obj.xmax = xmax;
            obj.ymin = ymin;
            obj.ymax = ymax;
            if ~isempty(obj.params.userBkgndROI)
                obj.params.userBkgndROI = obj.params.userBkgndROI((ymin - offsets(1)):(ymax - offsets(1)), ...
                    (xmin - offsets(2)):(xmax - offsets(2)),:,:);
            end
            
        end
           
        function obj = drawCircles(obj, circleSize, varargin)
            % Allow user to click multiple circle ROIs in a montage view
            % usage:
            % map = map.drawCircles(circleSize)
            % Displays the F0 of the activityMap object in montage view, allow
            % user to click several times to draw multiple circular ROIs,
            % save the circular ROIs as binary masks (1 for inside, 0 for outside) 
            % in the 'ROIs' field of the activityMap object
            % circleSize is the radius of the circle in microns (the
            % pixelCal field of the activityMap object is used to convert
            % to pixels)
            %
            % Optional parameters
            % map = map.drawCircles(circleSize, channel)
            % specify which channel of the F0 to draw
            % default is 1 (green)
            %
            % map = map.drawCircles(circleSize, channel, replaceExistingROIs)
            % if replaceExistingROIs==1, this will overwrite the existing
            % ROIs saved in the activityMap object
            % default is replaceExistingROIs = 0 (i.e., add ROIs to
            % existing ROIs)
            
            replaceExistingROIs = false;
            channel = 1;
            switch length(varargin)
                case 1
                    channel = varargin{1};
                case 2
                    channel = varargin{1};
                    replaceExistingROIs = varargin{2};
            end
            channelForDisplay = obj.f0unfiltered(:,:,:,channel);
            [montageForDisplay, montageDims] = stackToMontage(permute(channelForDisplay,[2 1 4 3]));

            % ask user to select the area
            fig = figure('Name', 'Click to draw circles (see console for instructions)');
            imagesc(montageForDisplay);
            axis equal, axis tight;
            if replaceExistingROIs
                circlesMontage = clickCircles(fig, gca, circleSize/obj.rowCal); % rowCal is microns per pixel
            else
                numROIs = size(obj.ROIs,4);
                montageROIs = zeros(size(montageForDisplay,1), size(montageForDisplay,2), numROIs);
                if isempty(obj.ROIs)
                    montageROIs = [];
                else
                    for i=1:numROIs
                        montageROIs(:,:,i) = stackToMontage(permute(obj.ROIs(:,:,:,i),[2 1 4 3]));
                    end
                end
                circlesMontage = clickCircles(fig, gca, circleSize/obj.rowCal, montageROIs);
            end
            result = montageToStack(circlesMontage,montageDims);
            % transpose image back into how Matlab stores the original data
            result = permute(result, [2 1 4 3]);
            result = logical(result);
            close(fig);
            obj.ROIs = result;

        end
        
        function [] = displayROIs(obj)
            figure;
            compressROIs = zeros(size(obj.ROIs,1),size(obj.ROIs,2),size(obj.ROIs,3));
            for i=1:size(obj.ROIs,1)
                for j=1:size(obj.ROIs,2)
                    for k=1:size(obj.ROIs,3)
                        
                        index = find(obj.ROIs(i,j,k,:));
                        if ~isempty(index)
%                             index
                            compressROIs(i,j,k) = index(1);
                            
                        end
                    end
                end
            end
            overlayImage = overlayColorOnGrayscale(permute(compressROIs, [2 1 3 4]), ...
                    permute(obj.f0(:,:,:,1), [2 1 3 4]), [0 max(compressROIs(:))], 0);
%             [montageForDisplay, montageDims] = stackToMontage(permute(compressROIs,[2 1 4 3]));
            imagesc(stackToMontage(overlayImage));
            axis equal, axis tight;
            
            
        end
        
        function obj = calcROITimeSeries(obj)
            % read the movie
            fprintf(strcat('Loading file ', obj.params.pathName, obj.params.fileName, '...\n'));
            strcat(obj.params.pathName,obj.params.fileName)
            [movie, ~] = readScanImageTiffLinLab(strcat(obj.params.pathName,obj.params.fileName));
            movie = single(movie);
 
            % have not written anything to deal with multiple channels.
            % figure this out later
            channel = 1;
            dffChannel = squeeze(movie(:,:,:,channel,:));
            dffChannel = dffChannel-obj.params.bkgnd(channel);
            % apply offsets in case the activityMap object has been
            % registered
            dffChannel = dffChannel((obj.ymin - obj.offsets(1)):(obj.ymax - obj.offsets(1)), ...
                (obj.xmin - obj.offsets(2)):(obj.xmax - obj.offsets(2)),:,:);

            % calculate dF/F for each ROI
            numROIs = size(obj.ROIs,4);
            obj.ROITimeSeries = zeros(numROIs, obj.nFrames);
            if size(obj.ROIs(:,:,:,1))~=size(dffChannel(:,:,:,1))
                error('size of ROI mask does not match size of movie');
            end
            
            ROImatrixsize = size(obj.ROIs);
            % reshape ROIs into a (x*y*z) x #ROIs matrix
            ROIslinear = reshape(obj.ROIs, prod(ROImatrixsize(1:3)), numROIs);
            % reshape the data into a (x*y*z) x timepoints matrix
            channelLinear = reshape(dffChannel, prod(ROImatrixsize(1:3)), size(dffChannel,4));
            
            obj.intraROIcorr = zeros(numROIs,1);
            for x=1:numROIs
                thisROIovertime = channelLinear(ROIslinear(:,x)==1,:);
                obj.ROITimeSeries(x,:) = mean(thisROIovertime,1);
                % correlate each pixel's time course with every other pixel
                % and get the mean of all intra-ROI inter-pixel
                % correlations
                if ~isempty(thisROIovertime)
                    obj.intraROIcorr(x) = meanExclDiag(corr(thisROIovertime'));
                    
                else
                    obj.intraROIcorr(x) = NaN;
                end
            end
%             for t=1:obj.nFrames
%                 thisFrame = dffChannel(:,:,:,t);
%                 for x=1:numROIs
%                     obj.ROITimeSeries(x,t) = mean(thisFrame(obj.ROIs(:,:,:,x)==1));
%                 end
%             end
            ROIF0 = mean(obj.ROITimeSeries(:,obj.prePeriod),2);
            
            for t = 1:obj.nFrames
                obj.ROITimeSeries(:,t) = (obj.ROITimeSeries(:,t) - ROIF0)./ROIF0;
            end
            
        end
        
        function obj = calcROIStimPeriod(obj)
            result = zeros(size(obj.ROIs,4), obj.numStimPeriods);
            if iscell(obj.stimPeriod)
                for i=1:obj.numStimPeriods
                    result(:,i) = mean(obj.ROITimeSeries(:,obj.stimPeriod{i}),2);
                end
            else
                result = mean(obj.ROITimeSeries(:,obj.stimPeriod),2);
            end
            obj.ROIStimPeriod = result;
                    
        end
        
        function result = ROIsignifResp(obj)
            stdPre = std(obj.ROITimeSeries(:,obj.prePeriod),1,2);
            result = (obj.ROIStimPeriod > obj.params.stdevThresMult*repmat(stdPre,1,obj.numStimPeriods));
        end

    end
    
    methods (Static)
        function result = corrOverTime(movie, f0, cap)
            % Returns the correlation between every individual frame of movie (an N-D
            % matrix, i.e. a movie) and f0 (an (N-1)-D matrix meant to be the average
            % over time of the movie over a certain period of time)
            % Caps pixel intensities at a certain threshold defined in the
            % parameters: params.capForCorrCalc (default is 10)
            % NB: this is static because the movie has to be passed in as 
            % a parameter, because the full movie is not saved in the object, 
            % to avoid wasting memory
            dim = size(movie);
            numDims = length(dim);
            numFrames = dim(numDims);
            
            result = zeros(numFrames,1);
            
            movie(movie>cap) = cap;
            f0(f0>cap) = cap;
            
            for i=1:numFrames
                index = repmat({':'},1,numDims);
                index{numDims} = i;
                thisFrame = movie(index{:}); % e.g. if numDims is 3, this would be movie(:,:,i)
                
                result(i) = corr(thisFrame(:), f0(:));
            end
        end
 
        function result = outlineObject(image, channel)
            % Allow user to manually outline the object in a montage view
            % usage:
            % mask = outlineObject(image, channel);
            % Displays the specified channel of image in montage view, allow
            % user to draw multiple ROIs to define a 3D object, return a
            % binary mask (1 for in the object, 0 for outside the object)
            % Note: image must be 4D (x, y, z, channel)
            % Example usage to draw mask based on F0 of channel 1 of an
            % activityMap object:
            % mask = activityMap.outlineObject(map.f0, 1);
            channelForDisplay = image(:,:,:,channel);
            [montageForDisplay, montageDims] = stackToMontage(permute(channelForDisplay,[2 1 4 3]));

            % ask user to select the area
            fig = figure('Name', 'Draw the 3D ROI (see console for instructions)');
            min(montageForDisplay(:))
            max(montageForDisplay(:))
            montageForDisplayMin=min(montageForDisplay(:));
            montageForDisplayMax=max(montageForDisplay(:));
            imagesc(montageForDisplay, [montageForDisplayMin montageForDisplayMax/5]);
            axis equal, axis tight;
            maskMontage = multiROI(fig, gca);
            
            result = montageToStack(maskMontage,montageDims);
            % transpose image back into how Matlab stores the original data
            result = permute(result, [2 1 3]);
            result = logical(result);
            close(fig);
        end
        
       
    end
    
end