# Astro2

This Git folder is intended to serve as the holder for the DRDC Ottawa Astro/Photometry Open Source Software project. 

### V1 - Created/ To Be Created
1. Automatic Plate Solving (Using Pinpoint)
    1. Astroquery to be integrated (To avoid use of pinpoint)
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

### Recently added functions/Capabilities (Jack): Sorted most recent first
