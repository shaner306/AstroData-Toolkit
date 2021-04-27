# Astro2

This Git folder is intended to serve as the holder for the DRDC Ottawa Astro/Photometry Open Source Software project. 
## Version 1.0
### AstroSolver (In Progress)
1. Stage 1 - Data Input 
    - Enter Image Folders (Decide on Parameters, different filtered images)
    - Execute Program
    - Read Folder and FITS files, save HEADER data of first file
    - Generate Output Logs
    - Import Reference Star Database data for searching
 
2. Stage 2 - Data Processing
    - Automatic Plate Solving (Using Pinpoint)
    - Astroquery to be integrated (To avoid use of pinpoint)

3. Stage 3 - Data Output
### AstroReducer (In Progress)


#### Upcoming
    - Documentation and User Manual
    - Comments
3. Basic Star Matching and Transform Calculations
    - With errors and Std dev integrated
4. Standard Magnitude Calculations of Basic Stars (To be completed)
5. Basic Star Data Collection and Output (Pinpoint Solving Point Stars)
      - FWHM
      - Magnitudes
      - Instrumental Mag
      - Flux, Stellar Flux, Visual Magnitude, gaussian data, errors
      - Location

5. TRM mode - _TRMtester.py_
  	1. Streak Detection Method
        - Simple Method (Weighted Centroids/ Length filtering)
        - Complex Method (Matched Filter/ Fast Fourier Transform)

  	2. Background Extraction
        - Sextractor
        - Polynomial Form Fit
        - Mean Background
        - Filter and Box Size TBD
        
     3. Creating Stars file
        - ImageStarsCache stars Creation
        - Pinpoint Interpretation
        - Pinpoint Solving
        
 
6.   Light Curves
        - Satellites
        - Point Stars
        - Streaked Stars (Possibly)


7. Satellite Matching, Tracking, and Photometry
8. Data Outputs
    - Include success statistics (e.g. number of images that were quality data points)
    - Errors and Margins
    - Calculation Constants
    - Astrometric Data
    - Photometric Data
    - Astrodynamic Data
    - 

### Recently added functions/Capabilities (Shane): Sorted most recent first
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
