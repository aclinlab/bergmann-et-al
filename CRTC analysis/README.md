# Explanation of CRTC analysis code

These are ImageJ macros for measuring nuclear localization of CRTC-GFP. To use, add them to the Startup Macros file for ImageJ. (Note there is a utility function `tagROIs`.)

## `CRTCAnalysis`

* Clear the ROI manager before starting
* Draw manual ROIs in triplets, and add each one in sequence to the ROI manager (each cell is a triplet).
	1. The whole cell
	2. The nucleus
	3. A background ROI
* Run `CRTCAnalysis`
* This will measure the green fluorescence inside the nucleus (ROI2), outside the nucleus (ROI1 XOR ROI2), and the background (ROI3) - and will do this for all the cells analysed 
* The macro automatically determines the green channel as the "nChannels-1"th channel
So if there are 3 channels, it takes the 2nd channel
If there are 2 channels, it takes the 1st channel
* It will output the results to the Log window like
`277.3569	215.6654	73.8795	6.9490	9.4530`
(That is, one row per cell)
* These numbers are in the order: 
`inMean	exMean	bkgndMean	inArea	exArea`
* The numbers will have been copied to the clipboard for pasting into a spreadsheet
* The macro will also auto-save the ROIs as a .zip file. If your file is myFile.nd2, it will save the ROIs as myFile-ROIs-1.zip, in the same directory as your original file. The -1 will change to -2, -3, etc., if there is already an ROIs zip file in the directory for the same image file. 

## `CRTCAnalysisOneBkgnd`

Same as above, except you don't draw a separate background ROI for every cell. The first ROI is the background, and the same ROI will be applied as the background to every cell (even across different slices - so the background ROI needs to be properly background across all z-slices).

Then draw the ROIs for each cell in pairs:
	1. The whole cell
	2. The nucleus
