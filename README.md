# AstroSolve

## Before you Start
1. Download Astrosolver repo - V1 AstroSolver folder
2. Download UCAC4 Catalog - http://forums.dc3.com/showthread.php?4694-Downloading-and-Using-PinPoint-Reference-Catalogs-(Updated-April-2021)
3. Have a registered copy of Pinpoint Astrometric Engine installed http://pinpoint.dc3.com
4. Have the properly formatted reference star .csv file - Reference Docs Folder


## GUI Workflow and Inputs
To run the Image Processor “AstroSolver”, run the GUI python script. The GUI will open a popup window displaying the following screen.


![AA354318-3DDD-4018-B6D3-0CE999D718E1](https://user-images.githubusercontent.com/75094714/145273674-39e4da25-4b04-4802-be28-cbdc1fc8f64a.png)


This screen allows you to adjust the parameters and inputs that will be processed. 

There are two tabs;
1. AstroSolver - Solves the Image using pinpoint, and then conducts photometric processing and data extraction.
2. Image Reduction - Conduct the pre-processing image subtraction and reduction, for Darks, Biases, and Flats.

### AstroSolver Program

The AstroSolver Program has 3 primary inputs;
1. The Image Folder - Where your star or satellite images are contained. 
2. The Catalog Folder - Where the star catalog that you wish to use for star matching is contained. We use UCAC4 as the test catalog.
3. The Reference Star file - A CSV, or excel file containing the reference star data for a series of stars accross a particular sky region. Used to calculate and compare the standard magnitudes produced in the processor.

#### Pinpoint Solve and Save Data
Below the folder inputs are two checkboxes;
1. Pinpoint Solve - Use Pinpoint Software to plate solve the images for more accurate astrometric data. 
	1. Requires 64/32 bit Pinpoint installed
	2. Valid Pinpoint License
	3. Can only be run on Windows
	4. Does not work with Track Rate Mode images
2. Save Data - Selecting this option generates a folder that will be populated with the graphs, charts, and tables generated through the processor for later use. *Default is Selected.*

### Analysis Parameters
On the GUI under Analysis Parameters, there are two radio buttons for the various program configurations. Only one option in each box can be selected and should be based on the source of the Images and the Mode in which they were captured.

#### Source Capture Mode
There are two mode options;
1. Star Stare Mode (SSM) or Sidereal Stare Mode, - The camera is directed toward a stationary direction and the object passes through the frame.
	 * Stars appear as points, Satellites appear as streaks
2. Track Rate Mode - The camera is slewed at the same rate in which the satellite moves accross the sky. 
	* Stars appear as streaks, the satellite appears as a point source.

#### Image Source
There are two choices for Image Source;
1. Space-Based - Images are taken from a satellite in orbit.
2. Ground-Based - Images are taken from a ground source such as optical observatory, amateur CCD, etc.
This selection is self evident but ensure the correct source is selected as the calculations are widely innaccurate if the incorrect method is applied.

### Pinpoint Solve Parameters
/These parameters only apply if the Pinpoint Solve option is selected/

There are 8 parameters that can be adjusted for solving using Pinpoint;
1. Maximum Magnitude - The Max visual magnitude Pinpoint will attempt to scan the image for. 
The higher the number, the lower the brightness. 
If the number is too high pinpoint will detect too many objects and the solution will take much longer. If the number is too low, pinpoint may not detect any and will not solve.
*Default is 13*
2. Sigma Above Mean
3. Max Solve Time
4. Max Match Residual
5. Catalog Expansion
*Default is 0.8*
6. Catalog - Which star catalog pinpoint is using to conduct the star matching. 
*Default is UCAC4*
7. Use SExtractor - Instead of using pinpoints normally background estimation and removal, it will use the SExtractor algorithm. This should be selected if the point sources are dim and a higher level of accuracy is necessary.
*Default is NO*
9. All Sky Solve - Conducts all sky scan instead of differential plate solving. Uses astrometry.net and takes a long duration. Only use in the most difficult of plate solving, and not with large imagesets.
*Default is NO*

Once the Parameters have been set, click “Solve” and monitor the progress within the python console.

Once started the GUI will freeze as it is not updated in real time, while the processing takes place.

## Image Reduction
The second tab is the Image Reduction Tab.

### Inputs
There is only one input for the Image Reduction program;
1. The Image Folder - The main folder with subfolders for Bias, Darks, Flats, and Lights.

#### Creating Master Files
There are three selections for Image Reduction;
1. Use Darks
2. Use Flats
3. Use Biases

Normally all three options are selected, and a master file will be produced for all three image sets and applied to the light frames. However if master files already exist or the light frames have already been adjusted for one or more of the reduction types, the option can be deselected and the program will run by skipping that reduction step.

It is usually necessary to do Image Reduction prior to image processing, however if your images are already reduced, further image reduction may cause distortion and reduce the quality of the image and the point sources youre are attempting to analyze. 

### Space-Based (Experimental Feature)
For space-based Imagery, such as NEOSSat, select the space-based option. NEOSSat can not do Bias and Flat reduction, only Dark Subtraction, and requires the satellite to focus on a particular target.


# Instructions
1. Read the description of the parameters and inputs for the program above.
2. Determine if your Images must be reduced prior to processing. If YES, click the Image Reduction Tab. If NO, continue to step 3.
3. Browse and select the correct folders for each input. Ensure the folder location URL is formatted similarly to the image above.
4. Select Save Data checkbox in order to save your data to the image folder.
5. If the images require plate solving; select Pinpoint Solve. If not, deselect the option.
6. Select the appropriate Analysis Parameters, If you are trying to image a star, the mode is likely SSM, for a satellite it is likely TRM.
7. If Pinpoint Solve was selected in Step 4, adjust the parameters as necessary. *Unless you require particular settings, use the defaults for the first run.*
8. Click Solve.
9. Return to Python and monitor the console for the image processing program to run.

<hr>
# Notes
* We use pinpoint which is a paid software for astrometric plate solving, however we will explore other options such as astrometry.net and astropy libraries such as astroquery.
* There are no good free solutions for astrometic plate solving, especially for a large imageset.

## Upcoming in V2
- Integrate errors, std dev
- TRM mode - _TRMtester.py_
  	1. Streak Detection Method
        - Simple Method (Weighted Centroids/ Length filtering)
        - Complex Method (Matched Filter/ Fast Fourier Transform)
  	2. Background Extraction
        - Sextractor
        - Polynomial Form Fit
        - Mean Background
        - Filter and Box Size
     3. Creating Stars file
        - ImageStarsCache stars Creation
        - Pinpoint Interpretation
        - Pinpoint Solving
- Light Curves
        - Satellites
        - Point Stars
        - Streaked Stars (Possibly)
        
- Satellite Matching, Tracking, and Photometry
- Complex Data Outputs
    - Include success statistics (e.g. number of images that were quality data points)
    - Errors and Margins
    - Calculation Constants
    - Astrometric Data
    - Photometric Data
    - Astrodynamic Data

### Recently added functions/Capabilities (Shane): Sorted most recent first
- Gaussian and Moffat Fitting of PSF
- Background Estimators (Mean, Mode, Median, SExtractor, Filtered, _Polyfit Coming soon_) - April 26th
- Star Removal Background Generator + Star Eraser
- Basic track rate mode algorithm (tested on neossat only)
- Saturated star finder
- Iterative Background Extraction
- Edge clipping
- Moment, Eccentricity, Centroid, and Compact Calculator for streaks
- Point Source Flux Extraction (non pinpoint)
- Dedicated Functions file
- FAST Transforms Calculator of ref star image
- Reference star Search and Detection
- Github created
- Complex Star image data collection algorithm and gaussian fit (Squid inspired)
### Recently added functions/Capabilities (Jack): Sorted most recent first
#### Star processor
- Apply ground-based transforms to a set of test stars
- Calculate ground-based transforms
- Calculate space-based transforms
- Match stars in a reference file to stars detected in an image
- Convert star positions from x,y to RA/dec and Alt/Az
- Calculate flux, instrumental magnitude, and instrumental magnitude standard deviations
- Calculate FWHM of all sources in an image
- Detect point sources in an image
- Calculate median and standard deviation background values
- Read a fits file
#### Light curve/satellite processor
- Interpolate to find the time-resolved colour indices of each satellite
- Plot seeing (FWHM) as a function of time
- Have a UI to input the locations of the desired satellites
- Match detected sources to satellite positions
- Calculate statistics (fwhm, magnitude, standard deviation of the magnitude) on the detected sources
- Detect point sources and their x,y locations
