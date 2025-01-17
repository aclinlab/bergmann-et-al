/*
 * // this loop measures average CRTC-GFP fluorescence in the intra- and extra-nuclear ROIs
// ROIs should be in TRIPLETS
// 1st is the whole cell
// 2nd is the nucleus
// 3rd is the background
 */

macro "CRTCAnalysis" {

	name=getTitle();
	directory = getDirectory("image"); //File.directory;
	if (roiManager("count") < 1) exit("Please load CRTC ROIs generated by CRTCRoiCreator");
	NumberOfROI= (roiManager("count"))/3;
	if (floor(NumberOfROI)!=NumberOfROI) exit("number of ROIs in ROI manager is not an even multiple of 3");

	run("Clear Results");
	print("\\Clear"); // clear the log
	//print(name);
	//print(directory);
	selectImage(name);
	Stack.getDimensions(width, height, channels, slices, frames);

	channel = channels-1; // green channel is the second to last channel because red is the last channel
	//print(channel);
	inMean = newArray(NumberOfROI);
	inArea = newArray(NumberOfROI);
	exMean = newArray(NumberOfROI);
	exArea = newArray(NumberOfROI);
	bkgndMean = newArray(NumberOfROI);
	InExNucRatio = newArray(NumberOfROI);
	pxlperROI = newArray(NumberOfROI);
	PercentNucPxl = newArray(NumberOfROI);
	for (i = 0; i < NumberOfROI; i++) {
		selectImage(name);
		
		// intranuclear
		roiManager("select", ((i*3)+1));
		Stack.setChannel(channel); // green channel
		run("Measure");
		inMean[i] = getResult("Mean");
		inArea[i] = getResult("Area");
		// extranuclear (cytoplasmic)
		roiManager("select", newArray((i*3),(i*3)+1));
		roiManager("XOR"); // get the XOR of the whole cell and the nucleus
		Stack.setChannel(channel); // green channel
		run("Measure");
		exMean[i] = getResult("Mean");
		exArea[i] = getResult("Area");
		// bkgnd
		roiManager("select", ((i*3)+2));
		Stack.setChannel(channel); // green channel
		run("Measure");
		bkgndMean[i] = getResult("Mean");
	
		// whole cell just for debugging
		roiManager("select", (i*3));
		Stack.setChannel(channel); // green channel
		run("Measure");
	}	

	// Print
	String.resetBuffer;
	for (i=0; i<NumberOfROI; i++) {
		toPrint = d2s(inMean[i],4)+"	"+d2s(exMean[i],4)+"	"+d2s(bkgndMean[i],4)+"	"+d2s(inArea[i],4)+"	"+d2s(exArea[i],4);
		print(toPrint);
		String.append(toPrint+"\n");
	}
	String.copy(String.buffer);
	selectWindow("Log");

	// save ROIs
	prevROIFileExists = true;
	index = 1;
	while (prevROIFileExists==true) {
		ROIfilename = tagROI(name,index);
		//print(ROIfilename);
		if (File.exists(directory+ROIfilename)) {
			index++;
		} else {
			roiManager("Save",directory+ROIfilename);
			prevROIFileExists = false;
		}
	}
}

function tagROI(filename, index) {
	suffixIndex = lastIndexOf(filename,"."); // find the last '.' in the file name to find the suffix
	if (suffixIndex != -1) {
		namePrefix = substring(filename, 0, suffixIndex); // everything up to the .
	} else {
		namePrefix = filename;
	}
	return namePrefix+"-ROIs-"+index+".zip";
	

}


/*
 * // this loop measures average CRTC-GFP fluorescence in the intra- and extra-nuclear ROIs
// ROIs should be in PAIRS
// 1st is the whole cell
// 2nd is the nucleus
// THE FIRST ROI IS THE BACKGROUND - IT IS APPLIED ACROSS ALL SLICES
 */
 
 macro "CRTCAnalysisOneBkgnd" {

	name=getTitle();
	directory = getDirectory("image"); //File.directory;
	if (roiManager("count") < 1) exit("Please load CRTC ROIs generated by CRTCRoiCreator");
	numROIs = (roiManager("count"));
	numCells = (numROIs-1)/2;
	if (floor(numCells)!=numCells) exit("number of ROIs in ROI manager is not an odd number!");

	run("Clear Results");
	print("\\Clear"); // clear the log
	//print(name);
	//print(directory);
	selectImage(name);
	Stack.getDimensions(width, height, channels, slices, frames);

	channel = channels-1; // green channel is the second to last channel because red is the last channel
	//print(channel);
	inMean = newArray(numCells);
	inArea = newArray(numCells);
	exMean = newArray(numCells);
	exArea = newArray(numCells);
	bkgndMean = newArray(numCells);
	InExNucRatio = newArray(numCells);
	pxlperROI = newArray(numCells);
	PercentNucPxl = newArray(numCells);
	for (i = 0; i < numCells; i++) {
		selectImage(name);
		// intranuclear
		roiManager("select", ((i*2)+2)); // the 0th ROI is the background
		Stack.setChannel(channel); // green channel
		run("Measure");
		inMean[i] = getResult("Mean");
		inArea[i] = getResult("Area");
		// extranuclear (cytoplasmic)
		roiManager("select", newArray((i*2+1),(i*2)+2));
		roiManager("XOR"); // get the XOR of the whole cell and the nucleus
		Stack.setChannel(channel); // green channel
		run("Measure");
		exMean[i] = getResult("Mean");
		exArea[i] = getResult("Area");
		// store the current slice number to apply to the background
		currentSlice = getSliceNumber();
		// bkgnd
		roiManager("select", 0); // background is the first ROI
		setSlice(currentSlice);
		Stack.setChannel(channel); // green channel
		run("Measure");
		bkgndMean[i] = getResult("Mean");
	
		// whole cell just for debugging
		roiManager("select", (i*2+1));
		Stack.setChannel(channel); // green channel
		run("Measure");
	}	

	// Print
	String.resetBuffer;
	for (i=0; i<numCells; i++) {
		toPrint = d2s(inMean[i],4)+"	"+d2s(exMean[i],4)+"	"+d2s(bkgndMean[i],4)+"	"+d2s(inArea[i],4)+"	"+d2s(exArea[i],4);
		print(toPrint);
		String.append(toPrint+"\n");
	}
	String.copy(String.buffer);
	selectWindow("Log");

	// save ROIs
	prevROIFileExists = true;
	index = 1;
	while (prevROIFileExists==true) {
		ROIfilename = tagROI(name,index);
		//print(ROIfilename);
		if (File.exists(directory+ROIfilename)) {
			index++;
		} else {
			roiManager("Save",directory+ROIfilename);
			prevROIFileExists = false;
		}
	}
}