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
3. The Reference Star file - A CSV, or excel file containing the reference 
   star data for a series of stars across a particular sky region. Used to 
   calculate and compare the standard magnitudes produced in the processor.

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
2. Track Rate Mode - The camera is slewed at the same rate in which the 
   satellite moves across the sky. 
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
9. All Sky Solve - Conducts all sky scan instead of differential plate 
   solving. Uses astrometry.net and takes a long duration. Only use in the most difficult of plate solving, and not with large image sets.
*Default is NO*

Once the Parameters have been set, click “Solve” and monitor the progress within the python console.

Once started the GUI will freeze as it is not updated in real time, while the processing takes place.

## Image Reduction
The second tab is the Image Reduction Tab.

### Inputs
There are Four initial inputs for the Image Reduction program;
1. The Image Folder - The folder of images to be calibrated
2. Bias Frames - The folder containing the bias frames. Note: this 
   selection is optional: without bias frames the program will default to 
   not dark scale the images.
3. Dark Frames - The folder containing the dark frames
4. Flat Frames - The folder containing the flat frames

#### Creating Master Files
There are four selections for Image Reduction under the directory tab;
1. Create Darks
2. Create Flats
3. Create Biases (Optional)
4. Mask for Outliers

If  master files do not exist, the first three options should be enabled.
If master files already exist, enable image correction from existing 
masters and point the directory to the master file collection.

To only create master files without reuducing images, only provide inputs for the Create Darks, Flats and Bias without an Image Directory.

It is usually advised to do Image Reduction prior to image processing, 
however if your images are already reduced, further image reduction may 
cause distortion and reduce the quality of the image and the point sources 
you are attempting to analyze. 

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

# Modules - Function Order
## A. Image Reduction
## B. Pinpoint Solve
## C. Ground-Based Photometric Processing (SSM)
## D. Space-Based Photometric Processing (SSM)
## E. Satellite Photometric Processing (TRM)
## F. (Coming soon) Satellite Astrometric Processing (TRM)
## G. Data Visualization


