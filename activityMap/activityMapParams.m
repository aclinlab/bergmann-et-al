classdef activityMapParams
    % activityMapParams: This class holds the parameters needed to construct an activityMap object
    % Open activityMapParams.m to see what the parameters are and which need to be set 
    properties
        % these need to be set before constructing an activityMap
        
        preLimits{mustBeNumeric} % define time limits of pre-stimulus period for calculating F0, in seconds, 1x2 matrix (e.g. [15 20])
        stimLimits % define time limits of stimulus period for calculating dF/F, in seconds, 1x2 matrix (e.g. [25 30])
        ctrlLimits{mustBeNumeric} %@double % define time limits (in seconds) of a control period, where you want to subtract dF/F from the real dF/F, 1x2 matrix (e.g. [25 30])
        
        % these have default values but you might want to set them
        
        corrThres{mustBeNumeric} = 0 % activityMap calculates the correlation of each frame with F0. If the object moves a lot, the correlation will fall. Set corrThres to exclude frames that fall below a certain correlation value. Default 0
        stdevThresMult{mustBeNumeric} = 2 % each pixel has a certain standard deviation during the F0 period. Its deltaF must exceed stdevThresMult x that stdev to count as "responding". Default 2
        bkgnd{mustBeNumeric} = 0 % background value to subtract from all pixels. Default 0. NB bkgnd should be a nChannels x 1 array if you care about having different backgrounds in each channel
        capForCorrCalc{mustBeNumeric} = 10 % for excluding motion, correlation of each frame with F0 may change not from motion but from genuine responses. set capForCorrCalc to prevent real responses from lowering correlation
        filter = fspecial('gaussian',5,2) % what filter to smooth the movie with in xy. Default fspecial('gaussian',5,2)
        userDrawBkgnd{mustBeNumericOrLogical} = false % if true, activityMap will prompt user to draw a freehand ROI of what pixels to use to calculate the bkgnd. If false, activityMap will use the property bkgnd 
        userDrawBkgndChannel = 1 % what channel to display for user to draw the background ROI
        timeSmoothingWindow = 1
        
        % these can be set by the activityMap constructor function
        
        fileName % file name (string). Can be set by activityMap constructor by user selecting file
        pathName % path name (string). Can be set by activityMap constructor by user selecting file
        userBkgndROI % an ROI of the background.. Can be set by activityMap constructor by user drawing the ROI, but you can copy the ROI from somewhere else to save time
        
        % optional parameters
        
        description % Optional. Use this if you want to generate many maps at once and annotate what each one is
    end
    
    properties (SetAccess = immutable)
        % update software version at each "release" to the lab
        softwareVersion = '1.3' %21/2/2018
    end
    
    methods
        function obj = activityMapParams
        end
    end
end