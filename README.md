# Astro2

This Git folder is intended to serve as the holder for the DRDC Ottawa Astro/Photometry Open Source Software project. 

### V1 - Created/ To Be Created
1. Integrated Pinpoint Solving
2. Basic Star Matching and Transform Calculations
    - With errors and Std dev integrated
3. Standard Magnitude Calculations of Basic Stars (To be completed)
4. Basic Star Data Collection and Output (Pinpoint Solving Point Stars)
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
