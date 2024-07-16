# activityMap

Type ```help activityMap``` or ```doc activityMap``` at the Matlab prompt to see the instructions for use.

Access the help files of the functions within the activityMap class using, e.g.:

```
help activityMap.overlayMontage
```

### Instructions and sample commands:

#### 1. Set parameters

```params = activityMapParams;```

This creates an empty activityMapParams object with default values. **NB you don't have to name your activityMapParams object ```params``` - you can call it whatever you want.**

               preLimits: []
              stimLimits: []
              ctrlLimits: []
               corrThres: 0
          stdevThresMult: 2
                   bkgnd: 0
          capForCorrCalc: 10
                  filter: [5×5 double]
           userDrawBkgnd: 0
    userDrawBkgndChannel: 1
                fileName: ''
                pathName: ''
            userBkgndROI: []
             description: ''
         softwareVersion: '1.3'
         
You MUST set preLimits and stimLimits before using the params object to construct an activityMap object. These are the time limits of the pre-stimulus period (to be used for calculating F0) and the stimulus period (to be used for calculating dF/F). They are in seconds and should each be a 1x2 matrix (eg ```[15 20]``` meaning the period from 15-20 s of the movie).

```
params.preLimits = [0 15];
params.stimLimits = [15 20];
```

#### 2. Create an activityMap object:

```map = activityMap(params);```

This calls the constructor function of activityMap, which will read in a movie and calculate an activity map. **Again, you don't have to name it ```map```, you can call it whatever you want.** If you didn't define the file and path of the .tif file in your activityMapParams object, Matlab will prompt you to select a file. After it's done calculating, Matlab will print something like this out:

    Loading file/Users/aclin/Google Drive/Lab/Data/Hoger/171023f1_00002.tif...
    Total time to read in file was: 
    Elapsed time is 4.900652 seconds.
    Filtering...
    Calculating DF...

    map = 

     activityMap with properties:

              f0: [4-D single]
         f0stdev: [4-D single]
             dff: [4-D single]
      dffNoThres: [4-D single]
    correlations: [225x1 double]
      mismatches: [0x1 double]
       prePeriod: [1x75 double]
      stimPeriod: [1x36 double]
      ctrlPeriod: []
      dimensions: [128 128 18 2 225]
       frameRate: 5.0005
           nRows: 128
        nColumns: 128
         nSlices: 18
       nChannels: 2
         nFrames: 225
          params: [1x1 activityMapParams]

#### 3. Visualize the map

To visualise the map as a montage where the dF/F (false-colored) is overlaid on the F0, at the prompt type something like:

```figure,imagesc(map.overlayMontage());```

Read the documentation of overlayMontage for more options. E.g. to give a fixed display range and overlay the dF/F of channel 1 on the F0 of channel 2, type:

```figure,imagesc(map.overlayMontage([0 8], [1 2]));```

To enforce square pixels, type:

```axis equal, axis tight```

#### 4. Save the map

```map.saveMap();```

If your movie file was fileName.tif, your activityMap object will be saved in a file called fileName-activitymap.mat, as a variable called 'obj'. The .mat file will go in the current directory

#### 5. Run skeleton analysis

```map = map.makeSkeleton(2,20);```

This uses channel 2 to define the skeleton with a node spacing of 20 microns (Type ```help activityMap.makeSkeleton``` and `help skeleton` for more details. `makeSkeleton` carries out the instructions specified by `help skeleton`).

If you want to re-use masks from a previous activity map, you can directly/separately use the commands that `makeSkeleton` calls:

```
mask = activityMap.outlineObject(obj.f0, channel);
obj.skel = skeleton(mask, obj.pixelCal);
obj.skel = obj.skel.defineStartingPoint(obj.f0(:,:,:,channel));
obj.skel = obj.skel.createSpacedNodes(spacing);
obj.skel = obj.skel.createVoronoiMask();
obj.skel.drawSkeleton();
```

```map = map.skeletonAnalysis(1);```

This uses channel 1 to calculate dF/F for each evenly spaced division in the object


##### 5.1 Rotating movies for display 3D objects
To save a rotating version of a skeleton drawing, FIRST create a Matlab figure with the skeleton drawing (e.g. `map.skel.drawSkeleton()` or `map.skel.drawSkeleton(map.voronoiMeanResp)`), and make sure it is the figure on top, then run:

`skeleton.saveRotatingMap(azimuthRange, azimuthStepSize, elevation, fileName);`

with values replacing the parameter variable names given here. Type `help skeleton.saveRotatingMap` for details.


### Version history

1.2:
* Added saveMap() function to activityMap to let you save the activityMap object as a .mat file
* Added plotCorrelations() function to activityMap to let you plot the correlations for each frame vs F0 - to look for motion artifacts that you can remove by setting the corrThres parameters appropriately

1.3
* Added skeletonisation
